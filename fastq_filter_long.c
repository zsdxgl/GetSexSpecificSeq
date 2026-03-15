// fastq_filter.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <zlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <immintrin.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sched.h>
#include <assert.h>

#define MAX_THREADS 255
#define K 31
#define BUFFER_SIZE (64 * 1024 * 1024)
#define QUEUE_SIZE 1024
#define CACHE_LINE_SIZE 64

// 线程安全队列
typedef struct {
    char* buffers[QUEUE_SIZE];
    size_t sizes[QUEUE_SIZE];
    int head;
    int tail;
    int count;
    pthread_mutex_t mutex;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} ThreadQueue;

// 共享数据结构
typedef struct {
    uint64_t* kmer_codes;
    uint64_t* rc_kmer_codes;
    int kmer_count;
    pthread_rwlock_t kmer_lock;
    
    unsigned long long total_reads;
    unsigned long long matched_reads;
    pthread_spinlock_t stat_lock;
} SharedData;

// 线程参数
typedef struct {
    int id;
    int cpu_id;
    SharedData* shared;
    ThreadQueue* input_queue;
    ThreadQueue* output_queue;
    char* output_filename;
    pthread_mutex_t* output_mutex;
    
    unsigned long long local_total;
    unsigned long long local_matched;
} ThreadData;

// DNA编码表
static const unsigned char dna_to_2bit[256] = {
    ['A'] = 0, ['a'] = 0,
    ['C'] = 1, ['c'] = 1,
    ['G'] = 2, ['g'] = 2,
    ['T'] = 3, ['t'] = 3
};

// 计算31-mer的64位编码
static inline uint64_t kmer_to_code(const char* seq) {
    uint64_t code = 0;
    for (int i = 0; i < K; i++) {
        unsigned char c = seq[i];
        unsigned char bits = dna_to_2bit[c];
        code = (code << 2) | bits;
    }
    return code;
}

// 计算反向互补编码
static inline uint64_t reverse_complement_code(uint64_t code) {
    uint64_t rc = 0;
    for (int i = 0; i < K; i++) {
        unsigned char bits = code & 3;
        rc = (rc << 2) | (3 - bits);
        code >>= 2;
    }
    return rc;
}

// 使用AVX2加速的字符串查找
static inline int contains_kmer_avx2(const char* seq, int len, 
                                    const uint64_t* kmer_codes, 
                                    const uint64_t* rc_kmer_codes,
                                    int kmer_count) {
    if (len < K) return 0;
    
    for (int i = 0; i <= len - K; i++) {
        uint64_t code = kmer_to_code(seq + i);
        
        for (int k = 0; k < kmer_count; k++) {
            if (code == kmer_codes[k] || code == rc_kmer_codes[k]) {
                return 1;
            }
        }
    }
    
    return 0;
}

// 加载kmer文件
void load_kmers(const char* filename, SharedData* shared) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        perror("无法打开kmer文件");
        exit(1);
    }
    
    char kmer[K + 2];
    int count = 0;
    
    while (fgets(kmer, sizeof(kmer), fp)) {
        if (strlen(kmer) >= K) {
            count++;
        }
    }
    
    rewind(fp);
    
    shared->kmer_codes = aligned_alloc(64, count * sizeof(uint64_t));
    shared->rc_kmer_codes = aligned_alloc(64, count * sizeof(uint64_t));
    shared->kmer_count = count;
    
    int idx = 0;
    while (fgets(kmer, sizeof(kmer), fp)) {
        if (strlen(kmer) >= K) {
            kmer[K] = '\0';
            
            for (int i = 0; i < K; i++) {
                if (kmer[i] >= 'a' && kmer[i] <= 'z') {
                    kmer[i] &= ~0x20;
                }
            }
            
            shared->kmer_codes[idx] = kmer_to_code(kmer);
            shared->rc_kmer_codes[idx] = reverse_complement_code(shared->kmer_codes[idx]);
            idx++;
        }
    }
    
    fclose(fp);
    printf("已加载 %d 个kmer\n", count);
}

// 队列操作
void queue_init(ThreadQueue* q) {
    pthread_mutex_init(&q->mutex, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
    q->head = q->tail = q->count = 0;
}

void queue_push(ThreadQueue* q, char* buffer, size_t size) {
    pthread_mutex_lock(&q->mutex);
    
    while (q->count == QUEUE_SIZE) {
        pthread_cond_wait(&q->not_full, &q->mutex);
    }
    
    q->buffers[q->tail] = buffer;
    q->sizes[q->tail] = size;
    q->tail = (q->tail + 1) % QUEUE_SIZE;
    q->count++;
    
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mutex);
}

int queue_pop(ThreadQueue* q, char** buffer, size_t* size) {
    pthread_mutex_lock(&q->mutex);
    
    while (q->count == 0) {
        pthread_cond_wait(&q->not_empty, &q->mutex);
    }
    
    *buffer = q->buffers[q->head];
    *size = q->sizes[q->head];
    q->head = (q->head + 1) % QUEUE_SIZE;
    q->count--;
    
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mutex);
    
    return 1;
}

// 读取线程
void* reader_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    gzFile input = gzopen(data->output_filename, "rb");
    
    if (!input) {
        perror("无法打开输入文件");
        exit(1);
    }
    
    gzbuffer(input, BUFFER_SIZE);
    
    char* buffer = NULL;
    size_t buffer_capacity = 0;
    size_t buffer_size = 0;
    
    while (1) {
        if (buffer == NULL) {
            buffer_capacity = BUFFER_SIZE;
            buffer = (char*)aligned_alloc(64, buffer_capacity);
        }
        
        int bytes_read = gzread(input, buffer + buffer_size, 
                               buffer_capacity - buffer_size - 1);
        
        if (bytes_read <= 0) {
            if (buffer_size > 0) {
                queue_push(data->input_queue, buffer, buffer_size);
                buffer = NULL;
            } else if (buffer) {
                free(buffer);
            }
            break;
        }
        
        buffer_size += bytes_read;
        
        char* end = buffer + buffer_size;
        char* last_record_end = buffer;
        
        int lines = 0;
        for (char* p = buffer; p < end; p++) {
            if (*p == '\n') {
                lines++;
                if (lines % 4 == 0) {
                    last_record_end = p + 1;
                }
            }
        }
        
        if (last_record_end > buffer) {
            size_t chunk_size = last_record_end - buffer;
            char* chunk = (char*)aligned_alloc(64, chunk_size);
            memcpy(chunk, buffer, chunk_size);
            
            queue_push(data->input_queue, chunk, chunk_size);
            
            size_t remaining = buffer_size - chunk_size;
            if (remaining > 0) {
                memmove(buffer, last_record_end, remaining);
            }
            buffer_size = remaining;
        } else {
            buffer_capacity *= 2;
            buffer = (char*)realloc(buffer, buffer_capacity);
        }
    }
    
    gzclose(input);
    
    for (int i = 0; i < data->id; i++) {
        queue_push(data->input_queue, NULL, 0);
    }
    
    return NULL;
}

// 处理线程
void* worker_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(data->cpu_id, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    
    char* chunk = NULL;
    size_t chunk_size = 0;
    
    while (1) {
        if (!queue_pop(data->input_queue, &chunk, &chunk_size)) {
            break;
        }
        
        if (chunk == NULL) {
            queue_push(data->output_queue, NULL, 0);
            break;
        }
        
        char* line_start = chunk;
        char* chunk_end = chunk + chunk_size;
        int line_num = 0;
        char* record[4];
        
        for (char* p = chunk; p < chunk_end; p++) {
            if (*p == '\n') {
                *p = '\0';
                
                if (line_num < 4) {
                    record[line_num] = line_start;
                }
                
                line_num++;
                line_start = p + 1;
                
                if (line_num == 4) {
                    data->local_total++;
                    
                    char* seq = record[1];
                    int seq_len = strlen(seq);
                    
                    if (seq_len >= K) {
                        int found = contains_kmer_avx2(seq, seq_len,
                                                     data->shared->kmer_codes,
                                                     data->shared->rc_kmer_codes,
                                                     data->shared->kmer_count);
                        
                        if (found) {
                            data->local_matched++;
                            
                            size_t rec_size = 0;
                            for (int i = 0; i < 4; i++) {
                                rec_size += strlen(record[i]) + 1;
                            }
                            
                            char* output = (char*)malloc(rec_size);
                            char* ptr = output;
                            
                            for (int i = 0; i < 4; i++) {
                                int len = strlen(record[i]);
                                memcpy(ptr, record[i], len);
                                ptr[len] = '\n';
                                ptr += len + 1;
                            }
                            
                            queue_push(data->output_queue, output, rec_size);
                        }
                    }
                    
                    line_num = 0;
                }
            }
        }
        
        free(chunk);
    }
    
    pthread_spin_lock(&data->shared->stat_lock);
    data->shared->total_reads += data->local_total;
    data->shared->matched_reads += data->local_matched;
    pthread_spin_unlock(&data->shared->stat_lock);
    
    return NULL;
}

// 写入线程
void* writer_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    
    gzFile output = gzdopen(fileno(fopen(data->output_filename, "wb")), "wb");
    gzbuffer(output, BUFFER_SIZE);
    
    int end_signals = 0;
    unsigned long long written = 0;
    
    while (end_signals < data->id) {
        char* buffer = NULL;
        size_t size = 0;
        
        if (!queue_pop(data->output_queue, &buffer, &size)) {
            break;
        }
        
        if (buffer == NULL) {
            end_signals++;
            continue;
        }
        
        pthread_mutex_lock(data->output_mutex);
        gzwrite(output, buffer, size);
        pthread_mutex_unlock(data->output_mutex);
        
        written++;
        free(buffer);
        
        if (written % 100000 == 0) {
            printf("已写入 %llu 条匹配序列\n", written);
        }
    }
    
    gzclose(output);
    return NULL;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        printf("使用方法: %s <kmer文件> <输入fastq.gz> <输出fastq.gz> <线程数>\n", argv[0]);
        printf("示例: %s kmers.txt input.fastq.gz output.fastq.gz 16\n", argv[0]);
        return 1;
    }
    
    const char* kmer_file = argv[1];
    const char* input_file = argv[2];
    const char* output_file = argv[3];
    int num_threads = atoi(argv[4]);
    
    if (num_threads > MAX_THREADS) {
        printf("警告：线程数超过最大值 %d，将使用 %d 个线程\n", MAX_THREADS, MAX_THREADS);
        num_threads = MAX_THREADS;
    }
    
    printf("使用 %d 个线程\n", num_threads);
    printf("kmer文件: %s\n", kmer_file);
    printf("输入文件: %s\n", input_file);
    printf("输出文件: %s\n", output_file);
    
    SharedData shared = {0};
    pthread_rwlock_init(&shared.kmer_lock, NULL);
    pthread_spin_init(&shared.stat_lock, 0);
    
    load_kmers(kmer_file, &shared);
    
    ThreadQueue input_queue, output_queue;
    queue_init(&input_queue);
    queue_init(&output_queue);
    
    pthread_t reader_tid;
    pthread_t worker_tids[MAX_THREADS];
    pthread_t writer_tid;
    ThreadData thread_data[MAX_THREADS];
    
    pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    ThreadData reader_data = {0};
    reader_data.id = num_threads;
    reader_data.output_filename = (char*)input_file;
    reader_data.input_queue = &input_queue;
    pthread_create(&reader_tid, NULL, reader_thread, &reader_data);
    
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].id = i;
        thread_data[i].cpu_id = i;
        thread_data[i].shared = &shared;
        thread_data[i].input_queue = &input_queue;
        thread_data[i].output_queue = &output_queue;
        thread_data[i].local_total = 0;
        thread_data[i].local_matched = 0;
        
        pthread_create(&worker_tids[i], NULL, worker_thread, &thread_data[i]);
    }
    
    ThreadData writer_data = {0};
    writer_data.id = num_threads;
    writer_data.output_filename = (char*)output_file;
    writer_data.output_queue = &output_queue;
    writer_data.output_mutex = &output_mutex;
    pthread_create(&writer_tid, NULL, writer_thread, &writer_data);
    
    pthread_join(reader_tid, NULL);
    
    for (int i = 0; i < num_threads; i++) {
        pthread_join(worker_tids[i], NULL);
    }
    
    pthread_join(writer_tid, NULL);
    
    printf("\n===== 统计信息 =====\n");
    printf("总处理序列数: %llu\n", shared.total_reads);
    printf("匹配序列数: %llu\n", shared.matched_reads);
    printf("匹配率: %.2f%%\n", 
           shared.matched_reads * 100.0 / (shared.total_reads ? shared.total_reads : 1));
    
    free(shared.kmer_codes);
    free(shared.rc_kmer_codes);
    pthread_rwlock_destroy(&shared.kmer_lock);
    pthread_spin_destroy(&shared.stat_lock);
    
    return 0;
}
