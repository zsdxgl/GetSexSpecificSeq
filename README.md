# GetSexSpecificSeq

1. get kmer
./getkmer.sh


less  getkmer.sh
./getkmer.pl 31
##output count.31.matrix

##
2.get sex specific kmer
awk '$3>=35 && $4==0 {print $0}' count.31.matrix >kmer_maleS.txt
cut -f 1 kmer_maleS.txt >kmer_maleS_seq.txt


###
3.get reads
##build
gcc -O3 -march=native -flto -o fastq_filter_short fastq_filter_short.c -lz -lpthread
gcc -O3 -march=native -flto -o fastq_filter_long fastq_filter_long.c -lz -lpthread
./paired_end_kmer_filter ./reseq/ kmerseq_maleS.txt maleS_1.fq.gz maleS_2.fq.gz 255
./getseq_bit -f reseq/ -k kmer_maleS_seq.txt -o kmer_maleS -t 36 
