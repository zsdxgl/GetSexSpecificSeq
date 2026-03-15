#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>

// 最大配置
#define MAX_THREADS 512
#define MAX_PATH_LEN 1024
#define MAX_SAMPLE_NAME 256
#define BUFFER_SIZE 262144
#define HASH_TABLE_SIZE 1000003
#define TASK_CHUNK_SIZE 10000000  // 每个任务块处理约1000万读数
#define MAX_QUEUE_SIZE 10000

// 原子操作
#define ATOMIC_ADD(ptr, val) __sync_add_and_fetch(ptr, val)
#define ATOMIC_INC(ptr) __sync_add_and_fetch(ptr, 1)

// 哈希表（快速查找）
typedef struct HashNode {
    uint64_t key;
    struct HashNode* next;
} HashNode;

typedef struct {
    HashNode** buckets;
    int size;
} HashTable;

// 样本信息
typedef struct {
    char r1_path[MAX_PATH_LEN];
    char r2_path[MAX_PATH_LEN];
    char sample_name[MAX_SAMPLE_NAME];
    uint64_t file_size_r1;
    uint64_t file_size_r2;
    uint64_t matches;
    uint64_t processed;
    uint64_t estimated_reads;
} SampleInfo;

// 任务块结构 - 将大文件分割成小块
typedef struct {
    int sample_idx;
    uint64_t start_byte_r1;
    uint64_t end_byte_r1;
    uint64_t start_byte_r2;
    uint64_t end_byte_r2;
    int chunk_id;
} TaskChunk;

// 线程安全的任务队列
typedef struct {
    TaskChunk* tasks;
    volatile int head;
    volatile int tail;
    int capacity;
    pthread_mutex_t lock;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} TaskQueue;

// 线程数据
typedef struct {
    int id;
    SampleInfo* samples;
    int total_samples;
    HashTable* kmer_hash;
    uint32_t kmer_len;
    uint64_t kmer_mask;
    gzFile out_r1;
    gzFile out_r2;
    pthread_mutex_t* output_lock;
    TaskQueue* task_queue;
    volatile uint64_t* thread_reads;
    volatile uint64_t* thread_matches;
    volatile int* should_exit;
} ThreadData;

// 全局统计
volatile uint64_t global_total_reads = 0;
volatile uint64_t global_total_matches = 0;
volatile uint64_t global_estimated_total_reads = 0;
volatile int global_samples_done = 0;
volatile int global_active_workers = 0;
volatile int global_total_chunks = 0;
volatile int global_chunks_done = 0;
pthread_mutex_t global_log_mutex = PTHREAD_MUTEX_INITIALIZER;
struct timeval global_start_time;

// 哈希表函数
HashTable* create_hash_table(int size) {
    HashTable* ht = malloc(sizeof(HashTable));
    ht->size = size;
    ht->buckets = calloc(size, sizeof(HashNode*));
    return ht;
}

void insert_kmer(HashTable* ht, uint64_t key) {
    int idx = key % ht->size;
    HashNode* node = malloc(sizeof(HashNode));
    node->key = key;
    node->next = ht->buckets[idx];
    ht->buckets[idx] = node;
}

bool contains_kmer(HashTable* ht, uint64_t key) {
    int idx = key % ht->size;
    HashNode* node = ht->buckets[idx];
    while (node) {
        if (node->key == key) return true;
        node = node->next;
    }
    return false;
}

void free_hash_table(HashTable* ht) {
    for (int i = 0; i < ht->size; i++) {
        HashNode* node = ht->buckets[i];
        while (node) {
            HashNode* next = node->next;
            free(node);
            node = next;
        }
    }
    free(ht->buckets);
    free(ht);
}

// 编码函数
int base_to_bits(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
}

uint64_t encode_kmer(const char* seq, int len) {
    uint64_t code = 0;
    for (int i = 0; i < len; i++) {
        code = (code << 2) | base_to_bits(seq[i]);
    }
    return code;
}

uint64_t reverse_complement_encode(uint64_t kmer, int len) {
    uint64_t rc = 0;
    for (int i = 0; i < len; i++) {
        int base = kmer & 3;
        int rc_base = 3 - base;
        rc = (rc << 2) | rc_base;
        kmer >>= 2;
    }
    return rc;
}

// 快速检查序列是否包含k-mer
bool fast_check_sequence(const char* seq, int seq_len, HashTable* ht, 
                        int k, uint64_t mask, uint64_t* matches) {
    if (seq_len < k) return false;
    
    int limit = seq_len - k;
    uint64_t code = encode_kmer(seq, k);
    
    for (int i = 0; i <= limit; i++) {
        if (i > 0) {
            int new_base = base_to_bits(seq[i + k - 1]);
            code = ((code << 2) | new_base) & mask;
        }
        
        if (contains_kmer(ht, code)) {
            if (matches) (*matches)++;
            return true;
        }
    }
    return false;
}

// 任务队列操作
TaskQueue* create_task_queue(int capacity) {
    TaskQueue* q = malloc(sizeof(TaskQueue));
    q->tasks = malloc(capacity * sizeof(TaskChunk));
    q->head = 0;
    q->tail = 0;
    q->capacity = capacity;
    pthread_mutex_init(&q->lock, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
    return q;
}

void push_task(TaskQueue* q, TaskChunk task) {
    pthread_mutex_lock(&q->lock);
    
    while ((q->tail + 1) % q->capacity == q->head) {
        pthread_cond_wait(&q->not_full, &q->lock);
    }
    
    q->tasks[q->tail] = task;
    q->tail = (q->tail + 1) % q->capacity;
    
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->lock);
}

TaskChunk pop_task(TaskQueue* q) {
    pthread_mutex_lock(&q->lock);
    
    while (q->head == q->tail) {
        pthread_cond_wait(&q->not_empty, &q->lock);
    }
    
    TaskChunk task = q->tasks[q->head];
    q->head = (q->head + 1) % q->capacity;
    
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->lock);
    
    return task;
}

int is_queue_empty(TaskQueue* q) {
    pthread_mutex_lock(&q->lock);
    int empty = (q->head == q->tail);
    pthread_mutex_unlock(&q->lock);
    return empty;
}

// 估计文件中的读数数量
uint64_t estimate_reads_from_file(const char* filename) {
    struct stat st;
    if (stat(filename, &st) == 0) {
        // 压缩FASTQ文件的典型大小：每个读数约300字节
        return st.st_size / 300;
    }
    return 1000000;  // 默认值
}

// 加载k-mer文件
int load_kmers_from_file(const char* filename, HashTable** ht_ptr, uint32_t* klen) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "错误: 无法打开k-mer文件 '%s'\n", filename);
        return -1;
    }
    
    char line[256];
    int count = 0;
    *klen = 0;
    
    // 第一遍：确定k-mer长度
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\r\n")] = 0;
        int len = strlen(line);
        if (len == 0) continue;
        if (*klen == 0) *klen = len;
        count++;
    }
    
    if (count == 0) {
        fclose(file);
        fprintf(stderr, "错误: k-mer文件为空\n");
        return -1;
    }
    
    rewind(file);
    HashTable* ht = create_hash_table(HASH_TABLE_SIZE);
    
    // 第二遍：加载k-mer
    int loaded = 0;
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\r\n")] = 0;
        if (strlen(line) != *klen) continue;
        
        uint64_t code = encode_kmer(line, *klen);
        uint64_t rc_code = reverse_complement_encode(code, *klen);
        insert_kmer(ht, code);
        insert_kmer(ht, rc_code);
        loaded++;
    }
    
    fclose(file);
    printf("已加载 %d 个k-mer（包含反向互补）\n", loaded * 2);
    
    *ht_ptr = ht;
    return loaded;
}

// 查找样本文件
int find_samples_in_dir(const char* dir_path, SampleInfo** samples_ptr) {
    DIR* dir = opendir(dir_path);
    if (!dir) {
        fprintf(stderr, "错误: 无法打开目录 '%s'\n", dir_path);
        return -1;
    }
    
    struct dirent* entry;
    SampleInfo* list = NULL;
    int capacity = 100;
    int count = 0;
    
    list = malloc(capacity * sizeof(SampleInfo));
    if (!list) {
        closedir(dir);
        return -1;
    }
    
    printf("扫描目录: %s\n", dir_path);
    while ((entry = readdir(dir)) != NULL) {
        char* fname = entry->d_name;
        if (!strstr(fname, ".fq.gz") && !strstr(fname, ".fastq.gz")) continue;
        
        // 检查R1文件
        char* suffix = NULL;
        if (strstr(fname, "_1.fq.gz")) suffix = "_1.fq.gz";
        else if (strstr(fname, "_R1.fq.gz")) suffix = "_R1.fq.gz";
        else if (strstr(fname, "_1.fastq.gz")) suffix = "_1.fastq.gz";
        else if (strstr(fname, "_R1.fastq.gz")) suffix = "_R1.fastq.gz";
        if (!suffix) continue;
        
        // 构建R2文件名
        char r2_name[MAX_PATH_LEN];
        strcpy(r2_name, fname);
        char* pos = strstr(r2_name, suffix);
        if (!pos) continue;
        
        if (strstr(suffix, "_1.")) strcpy(pos, "_2.fq.gz");
        else if (strstr(suffix, "_R1.")) strcpy(pos, "_R2.fq.gz");
        
        char r1_path[MAX_PATH_LEN], r2_path[MAX_PATH_LEN];
        snprintf(r1_path, sizeof(r1_path), "%s/%s", dir_path, fname);
        snprintf(r2_path, sizeof(r2_path), "%s/%s", dir_path, r2_name);
        
        if (access(r2_path, R_OK) != 0) {
            printf("警告: 缺少配对文件 %s\n", r2_name);
            continue;
        }
        
        // 提取样本名
        char sample_name[MAX_SAMPLE_NAME];
        strncpy(sample_name, fname, pos - r2_name);
        sample_name[pos - r2_name] = '\0';
        
        // 检查是否重复
        int exists = 0;
        for (int i = 0; i < count; i++) {
            if (strcmp(list[i].sample_name, sample_name) == 0) {
                exists = 1;
                break;
            }
        }
        if (exists) continue;
        
        // 添加到列表
        if (count >= capacity) {
            capacity *= 2;
            SampleInfo* new_list = realloc(list, capacity * sizeof(SampleInfo));
            if (!new_list) {
                free(list);
                closedir(dir);
                return -1;
            }
            list = new_list;
        }
        
        strcpy(list[count].r1_path, r1_path);
        strcpy(list[count].r2_path, r2_path);
        strcpy(list[count].sample_name, sample_name);
        
        // 获取文件大小
        struct stat st1, st2;
        if (stat(r1_path, &st1) == 0 && stat(r2_path, &st2) == 0) {
            list[count].file_size_r1 = st1.st_size;
            list[count].file_size_r2 = st2.st_size;
            list[count].estimated_reads = estimate_reads_from_file(r1_path);
        } else {
            list[count].file_size_r1 = 0;
            list[count].file_size_r2 = 0;
            list[count].estimated_reads = 1000000;
        }
        
        list[count].matches = 0;
        list[count].processed = 0;
        
        count++;
    }
    
    closedir(dir);
    
    if (count == 0) {
        free(list);
        fprintf(stderr, "错误: 未找到配对样本\n");
        return -1;
    }
    
    printf("找到 %d 个配对样本\n", count);
    *samples_ptr = list;
    return count;
}

// 为样本创建任务块
void create_tasks_for_sample(SampleInfo* sample, int sample_idx, 
                            TaskQueue* queue, uint64_t chunk_size_bytes) {
    
    uint64_t file_size = sample->file_size_r1;
    if (file_size == 0) {
        // 如果无法获取文件大小，创建一个默认任务
        TaskChunk task;
        task.sample_idx = sample_idx;
        task.chunk_id = 0;
        task.start_byte_r1 = 0;
        task.end_byte_r1 = UINT64_MAX;  // 处理整个文件
        task.start_byte_r2 = 0;
        task.end_byte_r2 = UINT64_MAX;
        push_task(queue, task);
        ATOMIC_INC(&global_total_chunks);
        return;
    }
    
    int num_chunks = (file_size + chunk_size_bytes - 1) / chunk_size_bytes;
    if (num_chunks == 0) num_chunks = 1;
    
    for (int i = 0; i < num_chunks; i++) {
        TaskChunk task;
        task.sample_idx = sample_idx;
        task.chunk_id = i;
        
        // 计算字节范围
        task.start_byte_r1 = i * chunk_size_bytes;
        task.end_byte_r1 = (i + 1) * chunk_size_bytes;
        if (task.end_byte_r1 > file_size) task.end_byte_r1 = file_size;
        
        // R2文件假设有相同的布局
        task.start_byte_r2 = task.start_byte_r1;
        task.end_byte_r2 = task.end_byte_r1;
        
        push_task(queue, task);
        ATOMIC_INC(&global_total_chunks);
    }
}

// 处理单个任务块
void process_task_chunk(TaskChunk task, SampleInfo* sample, HashTable* ht,
                       uint32_t k, uint64_t mask, gzFile out_r1, gzFile out_r2,
                       pthread_mutex_t* output_lock,
                       uint64_t* local_reads, uint64_t* local_matches) {
    
    gzFile in1 = gzopen(sample->r1_path, "rb");
    gzFile in2 = gzopen(sample->r2_path, "rb");
    
    if (!in1 || !in2) {
        fprintf(stderr, "错误: 无法打开文件块\n");
        return;
    }
    
    gzbuffer(in1, 1048576);
    gzbuffer(in2, 1048576);
    
    // 定位到块起始位置
    if (task.start_byte_r1 > 0) {
        if (gzseek(in1, task.start_byte_r1, SEEK_SET) == -1 ||
            gzseek(in2, task.start_byte_r2, SEEK_SET) == -1) {
            gzclose(in1);
            gzclose(in2);
            return;
        }
    }
    
    char buffer_r1[4][BUFFER_SIZE];
    char buffer_r2[4][BUFFER_SIZE];
    uint64_t bytes_processed = 0;
    uint64_t max_bytes = (task.end_byte_r1 == UINT64_MAX) ? UINT64_MAX : (task.end_byte_r1 - task.start_byte_r1);
    
    uint64_t temp_matches = 0;  // 将声明移到函数开头
    uint64_t read_count_since_last_update = 0;
    
    while (bytes_processed < max_bytes) {
        int lines_read = 0;
        
        // 读取一个完整的FASTQ记录
        for (int i = 0; i < 4; i++) {
            if (!gzgets(in1, buffer_r1[i], BUFFER_SIZE-1) ||
                !gzgets(in2, buffer_r2[i], BUFFER_SIZE-1)) {
                break;
            }
            buffer_r1[i][strcspn(buffer_r1[i], "\r\n")] = 0;
            buffer_r2[i][strcspn(buffer_r2[i], "\r\n")] = 0;
            bytes_processed += strlen(buffer_r1[i]) + 1;  // +1 for newline
            lines_read++;
        }
        
        if (lines_read != 4) break;
        
        (*local_reads)++;
        read_count_since_last_update++;
        
        // 检查是否包含k-mer
        int len1 = strlen(buffer_r1[1]);
        int len2 = strlen(buffer_r2[1]);
        bool match = false;
        
        if (len1 >= k) {
            match = fast_check_sequence(buffer_r1[1], len1, ht, k, mask, &temp_matches);
        }
        
        if (!match && len2 >= k) {
            match = fast_check_sequence(buffer_r2[1], len2, ht, k, mask, &temp_matches);
        }
        
        if (match) {
            (*local_matches)++;
            
            pthread_mutex_lock(output_lock);
            for (int i = 0; i < 4; i++) {
                gzprintf(out_r1, "%s\n", buffer_r1[i]);
                gzprintf(out_r2, "%s\n", buffer_r2[i]);
            }
            pthread_mutex_unlock(output_lock);
        }
        
        // 每处理10000个读数更新一次全局统计
        if (read_count_since_last_update >= 10000) {
            ATOMIC_ADD(&global_total_reads, read_count_since_last_update);
            ATOMIC_ADD(&global_total_matches, temp_matches);
            read_count_since_last_update = 0;
            temp_matches = 0;
        }
    }
    
    // 更新剩余统计
    if (read_count_since_last_update > 0) {
        ATOMIC_ADD(&global_total_reads, read_count_since_last_update);
        ATOMIC_ADD(&global_total_matches, temp_matches);
    }
    
    gzclose(in1);
    gzclose(in2);
}

// 工作线程函数
void* worker_thread_func(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    ATOMIC_INC(&global_active_workers);
    
    uint64_t local_reads = 0, local_matches = 0;
    
    while (!(*data->should_exit)) {
        // 尝试获取任务
        TaskChunk task;
        int got_task = 0;
        
        pthread_mutex_lock(&data->task_queue->lock);
        if (data->task_queue->head != data->task_queue->tail) {
            task = data->task_queue->tasks[data->task_queue->head];
            data->task_queue->head = (data->task_queue->head + 1) % data->task_queue->capacity;
            got_task = 1;
        }
        pthread_mutex_unlock(&data->task_queue->lock);
        
        if (!got_task) {
            // 队列为空，检查是否所有任务都已完成
            if (global_chunks_done >= global_total_chunks && global_total_chunks > 0) {
                break;
            }
            usleep(1000);  // 短暂等待
            continue;
        }
        
        // 处理任务块
        SampleInfo* sample = &data->samples[task.sample_idx];
        
        process_task_chunk(task, sample, data->kmer_hash, data->kmer_len,
                          data->kmer_mask, data->out_r1, data->out_r2,
                          data->output_lock, &local_reads, &local_matches);
        
        ATOMIC_INC(&global_chunks_done);
        
        // 更新线程本地统计
        data->thread_reads[data->id] += local_reads;
        data->thread_matches[data->id] += local_matches;
        local_reads = 0;
        local_matches = 0;
    }
    
    ATOMIC_ADD(&global_active_workers, -1);
    return NULL;
}

// 进度监控
void* progress_monitor_func(void* arg) {
    int* params = (int*)arg;
    int num_threads = params[0];
    int sample_count = params[1];
    
    struct timeval last_report, now;
    uint64_t last_reads = 0;
    int last_chunks = 0;
    
    gettimeofday(&last_report, NULL);
    sleep(5);  // 等待线程启动
    
    while (global_active_workers > 0 || global_chunks_done < global_total_chunks) {
        sleep(5);
        gettimeofday(&now, NULL);
        
        uint64_t reads = global_total_reads;
        uint64_t matches = global_total_matches;
        int chunks_done = global_chunks_done;
        int chunks_total = global_total_chunks;
        int active_workers = global_active_workers;
        
        double elapsed = (now.tv_sec - global_start_time.tv_sec) +
                        (now.tv_usec - global_start_time.tv_usec) / 1000000.0;
        
        double recent_speed = 0;
        if (elapsed > 5) {
            recent_speed = (reads - last_reads) / 5.0;
        }
        
        pthread_mutex_lock(&global_log_mutex);
        printf("\n=== 进度报告 [运行: %.0f秒] ===\n", elapsed);
        
        // 任务块进度
        if (chunks_total > 0) {
            printf("任务块: %d/%d (%.1f%%) 完成\n", 
                   chunks_done, chunks_total, (chunks_done*100.0)/chunks_total);
        }
        
        // 读数进度
        if (global_estimated_total_reads > 0) {
            double percent = (reads * 100.0) / global_estimated_total_reads;
            printf("读数处理: %lu/%lu (%.1f%%) 完成\n", 
                   reads, global_estimated_total_reads, percent);
        } else {
            printf("读数处理: %lu 个\n", reads);
        }
        
        // 匹配统计
        if (reads > 0) {
            printf("匹配: %lu (%.4f%%)\n", matches, (matches*100.0)/reads);
        }
        
        // 速度
        if (recent_speed > 0) {
            double avg_speed = reads / elapsed;
            printf("速度: 当前 %.0f 读数/秒, 平均 %.0f 读数/秒\n", 
                   recent_speed, avg_speed);
        }
        
        // 活跃线程
        printf("活跃线程: %d/%d\n", active_workers, num_threads);
        
        // 准确估算剩余时间
        if (recent_speed > 0 && global_estimated_total_reads > reads) {
            double remaining = global_estimated_total_reads - reads;
            double remaining_seconds = remaining / recent_speed;
            
            if (remaining_seconds < 60) {
                printf("预计剩余: %.0f 秒\n", remaining_seconds);
            } else if (remaining_seconds < 3600) {
                printf("预计剩余: %.1f 分钟\n", remaining_seconds / 60);
            } else {
                printf("预计剩余: %.1f 小时\n", remaining_seconds / 3600);
            }
        }
        
        printf("==================================\n");
        pthread_mutex_unlock(&global_log_mutex);
        
        last_reads = reads;
        last_chunks = chunks_done;
        last_report = now;
    }
    
    return NULL;
}

// 打印用法
void print_usage(char* prog_name) {
    printf("用法: %s [选项] <输入目录> <k-mer文件> <输出前缀>\n", prog_name);
    printf("选项:\n");
    printf("  -t <线程数>   指定使用的CPU核心数 (默认: 自动检测)\n");
    printf("  -c <块大小>   指定任务块大小(MB) (默认: 100, 即100MB)\n");
    printf("  -h           显示帮助信息\n");
    printf("\n示例:\n");
    printf("  %s -t 128 -c 200 ./data kmers.txt output\n", prog_name);
    printf("  %s ./data kmers.txt output\n", prog_name);
}

int main(int argc, char* argv[]) {
    printf("========================================\n");
    printf("fastq_filter_short - 极速k-mer序列提取工具\n");
    printf("工作窃取线程池 + 文件分块处理\n");
    printf("========================================\n\n");
    
    // 解析命令行参数
    char* seq_dir = NULL;
    char* kmer_file = NULL;
    char* output_prefix = NULL;
    int num_threads = sysconf(_SC_NPROCESSORS_ONLN);
    int chunk_size_mb = 100;  // 默认100MB每个块
    
    int opt;
    while ((opt = getopt(argc, argv, "t:c:h")) != -1) {
        switch (opt) {
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'c':
                chunk_size_mb = atoi(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                fprintf(stderr, "错误: 无效参数\n");
                print_usage(argv[0]);
                return 1;
        }
    }
    
    // 检查必需参数
    if (optind + 3 > argc) {
        fprintf(stderr, "错误: 缺少必需参数\n");
        print_usage(argv[0]);
        return 1;
    }
    
    seq_dir = argv[optind];
    kmer_file = argv[optind + 1];
    output_prefix = argv[optind + 2];
    
    // 限制最大线程数
    if (num_threads > MAX_THREADS) {
        printf("警告: 线程数限制为 %d\n", MAX_THREADS);
        num_threads = MAX_THREADS;
    }
    
    if (num_threads < 1) {
        num_threads = 1;
    }
    
    printf("配置:\n");
    printf("  输入目录: %s\n", seq_dir);
    printf("  k-mer文件: %s\n", kmer_file);
    printf("  输出前缀: %s\n", output_prefix);
    printf("  使用线程: %d\n", num_threads);
    printf("  任务块大小: %d MB\n", chunk_size_mb);
    printf("\n");
    
    gettimeofday(&global_start_time, NULL);
    
    // 1. 加载k-mer
    printf("步骤1: 加载k-mer...\n");
    HashTable* kmer_hash = NULL;
    uint32_t kmer_length = 0;
    
    if (load_kmers_from_file(kmer_file, &kmer_hash, &kmer_length) <= 0) {
        return 1;
    }
    
    uint64_t kmer_mask = (1ULL << (kmer_length * 2)) - 1;
    printf("k-mer长度: %u, 掩码: 0x%016lx\n", kmer_length, kmer_mask);
    
    // 2. 查找样本
    printf("步骤2: 扫描样本目录...\n");
    SampleInfo* samples = NULL;
    int sample_count = find_samples_in_dir(seq_dir, &samples);
    if (sample_count <= 0) {
        free_hash_table(kmer_hash);
        return 1;
    }
    
    // 估算总读数
    global_estimated_total_reads = 0;
    for (int i = 0; i < sample_count; i++) {
        global_estimated_total_reads += samples[i].estimated_reads;
    }
    printf("预估总读数: %lu\n", global_estimated_total_reads);
    
    // 3. 准备输出文件
    printf("步骤3: 准备输出文件...\n");
    char out_r1_name[MAX_PATH_LEN], out_r2_name[MAX_PATH_LEN];
    snprintf(out_r1_name, sizeof(out_r1_name), "%s_1.fq.gz", output_prefix);
    snprintf(out_r2_name, sizeof(out_r2_name), "%s_2.fq.gz", output_prefix);
    
    gzFile out_r1 = gzopen(out_r1_name, "wb");
    gzFile out_r2 = gzopen(out_r2_name, "wb");
    if (!out_r1 || !out_r2) {
        fprintf(stderr, "错误: 无法创建输出文件\n");
        free_hash_table(kmer_hash);
        free(samples);
        return 1;
    }
    
    gzbuffer(out_r1, 1048576);
    gzbuffer(out_r2, 1048576);
    pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    // 4. 创建任务队列并分割任务
    printf("步骤4: 创建任务队列...\n");
    TaskQueue* task_queue = create_task_queue(MAX_QUEUE_SIZE);
    uint64_t chunk_size_bytes = (uint64_t)chunk_size_mb * 1024LL * 1024LL;
    
    for (int i = 0; i < sample_count; i++) {
        create_tasks_for_sample(&samples[i], i, task_queue, chunk_size_bytes);
    }
    
    printf("创建了 %d 个任务块\n", global_total_chunks);
    
    // 5. 创建工作线程
    printf("步骤5: 启动 %d 个工作线程...\n", num_threads);
    pthread_t* threads = malloc(num_threads * sizeof(pthread_t));
    ThreadData* thread_data = malloc(num_threads * sizeof(ThreadData));
    pthread_t monitor_thread;
    
    volatile int should_exit = 0;
    volatile uint64_t* thread_reads = calloc(num_threads, sizeof(uint64_t));
    volatile uint64_t* thread_matches = calloc(num_threads, sizeof(uint64_t));
    
    int monitor_params[2] = {num_threads, sample_count};
    
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].id = i;
        thread_data[i].samples = samples;
        thread_data[i].total_samples = sample_count;
        thread_data[i].kmer_hash = kmer_hash;
        thread_data[i].kmer_len = kmer_length;
        thread_data[i].kmer_mask = kmer_mask;
        thread_data[i].out_r1 = out_r1;
        thread_data[i].out_r2 = out_r2;
        thread_data[i].output_lock = &output_mutex;
        thread_data[i].task_queue = task_queue;
        thread_data[i].thread_reads = thread_reads;
        thread_data[i].thread_matches = thread_matches;
        thread_data[i].should_exit = &should_exit;
        
        if (pthread_create(&threads[i], NULL, worker_thread_func, &thread_data[i]) != 0) {
            fprintf(stderr, "错误: 无法创建线程 %d\n", i);
            return 1;
        }
    }
    
    // 6. 启动进度监控
    pthread_create(&monitor_thread, NULL, progress_monitor_func, monitor_params);
    
    // 7. 等待所有线程完成
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    should_exit = 1;
    pthread_join(monitor_thread, NULL);
    
    // 8. 收尾工作
    gzclose(out_r1);
    gzclose(out_r2);
    pthread_mutex_destroy(&output_mutex);
    
    struct timeval end_time;
    gettimeofday(&end_time, NULL);
    double total_time = (end_time.tv_sec - global_start_time.tv_sec) +
                       (end_time.tv_usec - global_start_time.tv_usec) / 1000000.0;
    
    // 9. 最终报告
    printf("\n========================================\n");
    printf("处理完成！\n");
    printf("========================================\n");
    printf("样本数: %d\n", sample_count);
    printf("任务块数: %d\n", global_total_chunks);
    printf("总处理读数: %lu\n", global_total_reads);
    printf("总匹配读数: %lu (%.6f%%)\n", global_total_matches,
           global_total_reads ? (global_total_matches*100.0)/global_total_reads : 0.0);
    printf("总运行时间: %.2f 秒 (%.1f 分钟)\n", total_time, total_time/60.0);
    if (total_time > 0) {
        printf("平均速度: %.0f 读数/秒\n", global_total_reads / total_time);
    }
    printf("输出文件:\n  %s\n  %s\n", out_r1_name, out_r2_name);
    
    // 线程统计
    printf("\n线程统计:\n");
    uint64_t total_thread_reads = 0;
    for (int i = 0; i < num_threads; i++) {
        if (thread_reads[i] > 0) {
            printf("  线程 %d: %lu 读数, %lu 匹配\n", 
                   i, thread_reads[i], thread_matches[i]);
            total_thread_reads += thread_reads[i];
        }
    }
    
    // 样本统计
    printf("\n样本统计:\n");
    for (int i = 0; i < sample_count; i++) {
        if (samples[i].processed > 0) {
            printf("  样本 %s: %lu 读数, %lu 匹配 (%.4f%%)\n",
                   samples[i].sample_name, samples[i].processed, samples[i].matches,
                   (samples[i].matches * 100.0) / samples[i].processed);
        }
    }
    
    printf("========================================\n");
    
    // 10. 清理
    free_hash_table(kmer_hash);
    free(samples);
    free(threads);
    free(thread_data);
    free((void*)thread_reads);
    free((void*)thread_matches);
    free(task_queue->tasks);
    free(task_queue);
    
    return 0;
}
