# GetSexSpecificSeq

1. #get kmer
2. ./getkmer.sh
3. #less  getkmer.sh
4. #./getkmer.pl 31
5. ##output count.31.matrix
6. #get sex specific kmer
7. awk '$3>=35 && $4==0 {print $0}' count.31.matrix >kmer_maleS.txt
8. cut -f 1 kmer_maleS.txt >kmer_maleS_seq.txt
9. #get reads
10. ##build
11. gcc -O3 -march=native -flto -o fastq_filter_short fastq_filter_short.c -lz -lpthread
12. gcc -O3 -march=native -flto -o fastq_filter_long fastq_filter_long.c -lz -lpthread
13. ./paired_end_kmer_filter ./reseq/ kmerseq_maleS.txt maleS_1.fq.gz maleS_2.fq.gz 255
14. ./getseq_bit -f reseq/ -k kmer_maleS_seq.txt -o kmer_maleS -t 36 
