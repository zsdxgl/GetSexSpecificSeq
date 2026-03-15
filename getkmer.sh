#!/bin/sh
#SBATCH -J K31 #作业名，须更改
#SBATCH -N 1 #指定⼏个节点，此处规定设置为1
#SBATCH -c 10 #指定CPUs，此处规定不大于10
#SBATCH -p normal #指定Partition，选择normal或者inspur（默认）
#SBATCH -o %j.out #标准输出
#SBATCH -e %j.err #错误输出
#SBATCH -w node01 #可选node01、node02和node03
#conda activate kmtricks_env
./getkmer.pl 31

awk '$3>=35 && $4==0 {print $0}' count.31.matrix >kmer_maleS.txt

 makeblastdb -dbtype nucl -in ont_flye/assembly.fasta -input_type fasta -parse_seqids -out ont_flye
blastn -query kmer_maleS_seq.fa -db ont_flye  -out kmer_maleS.TO.ont_flye.txt -outfmt 6
sort -k9,9n kmer_maleS.TO.ont_flye.txt |awk '$7==1 && $8==31 && $3==100 {print}' |less -SN

cut -f 1 kmer_maleS.txt >kmer_maleS_seq.txt

##seqkit 太慢了，用fastq_filter


#利用C提取hifi序列，很快
#程序 fastq_filter.c
#编译 build.sh
运行 ./fastq_filter kmer_maleS_seq.txt 200 
nohup  ./fastq_filter kmer_maleS_seq.txt 200 &>fastq_filter.log &
echo 3277925 >pid.fastq_filter
#hifiasm  --b-cov 0 --h-cov -1 matching_reads.fastq.gz -l 0 -n 0 --ctg-n 0
hifiasm -o y_region -t 16 --hom-cov 28 --b-cov 0 --h-cov -1 -l 0 -n 7 --ctg-n 14 -u 1 --n-hap 1 matching_reads.fastq.gz
awk '/^S/{print ">"$2;print $3}' y_region.bp.p_ctg.gfa  > y_region.bp.p_ctg.fa
##illumina
#getset.c
#gcc -o getseq getseq.c -lz -lpthread -O3 -march=native

nohup ./getseq_bit -f reseq/ -k kmer_maleS_seq.txt -o kmer_maleS -t 36 &>getseq_bit.log &
echo 3278340 >pid.getseq_bit


##双末端序列提取,已经舍
gcc -O3 -march=native -flto -o fastq_filter_short fastq_filter_short.c -lz -lpthread
#gcc -O3 -march=native -pthread -lz paired_end_kmer_filter.c -o paired_end_kmer_filter
#./paired_end_kmer_filter ./reseq/ kmerseq_maleS.txt maleS_1.fq.gz maleS_2.fq.gz 255

####polish 

awk '$12>50 && ($4-$3)>0.8*$2 {print}' hifiraw2ont_fly.paf |cut -f 1 >hifiraw2ont_fly_selectedhifi.txt
seqtk subseq  hifi_kmer.fastq.gz hifiraw2ont_fly_selectedhifi.txt  |gzip >hifi_selected.fq.gz
minimap2 -x map-hifi ont_flye/assembly.fasta hifi_selected.fq.gz -a |samtools view -h -q 50  | samtools sort -O bam -o hifi_selected2ont_flye.bam
samtools index hifi_selected2ont_flye.bam
samtools consensus -f fastq -o ont_flye_assembly_consensus.fq.gz  -a -T ./ont_flye/assembly.fasta -@ 10 hifi_selected2ont_flye.bam
seqkit fq2fa ont_flye_assembly_consensus.fq.gz >ont_flye_assembly_consensus.fa

###valiate
minimap2 -x map-ont ont_flye_assembly_consensus.fa ont_kmer.fastq.gz >ontraw2ont_flye_assembly_consensus.paf
awk '$12>50 && $2*0.8<($4-$3) {print }' ontraw2ont_flye_assembly_consensus.paf |sort -k8,8n |less -SN
minimap2 -x map-hifi ont_flye_assembly_consensus.fa hifi_selected.fq.gz >hifi_selected2ont_flye_assembly_consensus.paf
awk '$12>50 && $2*0.8<($4-$3) {print }' hifi_selected2ont_flye_assembly_consensus.paf |sort -k8,8n |less -SN
##有效范围154975-1767156
makeblastdb -dbtype nucl -in ont_flye_assembly_consensus.fa -input_type fasta -parse_seqids -out ont_flye_assembly_consensus
blastn -query kmer_maleS_seq.fa -db ont_flye_assembly_consensus  -out kmer_maleS.TO.ont_flye_assembly_consensus.txt -outfmt 6
awk '$7==1 && $8==31 && $3==100 {print }' kmer_maleS.TO.ont_flye_assembly_consensus.txt |sort -k9,9n |less -SN
###有效范围173916-1747761
###replace the T2T.male.fa
minimap2 -x asm5 ../0.genome/T2T_male.fasta ont_flye_assembly_consensus.fa >ont_flye_assembly_consensus2T2T_male.paf
 cp ont_flye_assembly_consensus2T2T_male.paf replace_genome.paf
 perl ./replace_genome.pl  -r ../0.genome/T2T_male.fasta -q ont_flye_assembly_consensus.fa -a replace_genome.paf -o T2T_man_update.fa
minimap2 -x asm5 T2T_man_update.fa ont_flye_assembly_consensus.fa >ont_flye_assembly_consensus2T2T_man_update.paf

 
 makeblastdb -dbtype nucl -in T2T_man_update.fa -input_type fasta -parse_seqids -out T2T_man_update
blastn -query kmer_maleS_seq.fa -db T2T_man_update  -out kmer_maleS.TO.T2T_man_update.txt -outfmt 6
sort -k9,9n kmer_maleS.TO.T2T_man_update.txt |awk '$7==1 && $8==31 && $3==100 {print}' |less -SN
sort -k9,9n kmer_maleS.TO.T2T_man_update.txt |awk '$7==1 && $8==31 && $3==100 {print}' |awk '{print $2"\t"$9"\t"$2":"$9"\t"$11}' >kmer_maleS.gwas
minimap2 -x map-ont T2T_man_update.fa ont_kmer.fastq.gz >ontraw2T2T_man_update.paf
minimap2 -x map-hifi -t 200 T2T_man_update.fa hifi_selected.fq.gz >hifiraw2T2T_man_update.paf
awk '$12>50 && $2*0.8<($4-$3) {print }' hifi_selected2T2T_man_update.paf sort -k8,8n |less -N
minimap2 -x map-ont T2T_feman.fa ont_kmer.fastq.gz >ontraw2T2T_feman.paf
minimap2 -x map-hifi -t 200 T2T_female.fasta hifi_selected.fq.gz >hifiraw2T2T_female.paf

samtools faidx T2T_female.fasta chr10:35072671-37214343 >T2T_female.chr10-35072671-37214343.fa
blastn -query T2T_female.chr10-35072671-37214343.fa -db ont_flye_assembly_consensus  -out T2T_female.chr10-35072671-37214343.TO.ont_flye_assembly_consensus.txt -outfmt 6

cat T2T_female.chr10-35072671-37214343.fa   >female.male.seq.RC.fa 
RC.pl ont_flye_assembly_consensus.fa  >>female.male.seq.RC.fa
mafft --6merpair --thread -100 female.male.seq.RC.fa > female.male.seq.RC.mafft.fa


minimap2 -x asm5 -a T2T_man_update.fa T2T_female.chr10-35072671-37214343.  > T2Tfemaletarget2T2T_man_update.sam
samtools view -bS  T2Tfemaletarget2T2T_man_update.sam | samtools sort -o  T2Tfemaletarget2T2T_man_update.sorted.bam
samtools index  T2Tfemaletarget2T2T_man_update.sorted.bam
rm T2Tfemaletarget2T2T_man_update.sam
