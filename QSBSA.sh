#!/bin/bash
#SBATCH --job-name=BSA
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --mail-type=END
#SBATCH --mem=150G
#SBATCH --mail-user=qiaoshan.lin@uconn.edu
#SBATCH -o BSA_%j.out
#SBATCH -e BSA_%j.err

module load Trimmomatic/0.36
module load bowtie2/2.3.1 
module load samtools/1.3.1
module load bamtools/2.4.1
module load bcftools/1.3.1
module load R/3.4.3

#path to genome reference fasta file 
ref=
#path to forward reads fastq(.gz) file
read1=
#path to reverse reads fastq(.gz) file
read2=

mkdir tmp
cd tmp
####### reads trimming #######
echo Begin reads trimming...
java -jar $Trimmomatic PE \
 -phred33 -threads 8 \
 $read1 $read2 \
 FP.fq.gz FU.fq.gz \
 RP.fq.gz RU.fq.gz \
 LEADING:20 \
 TRAILING:20 \
 SLIDINGWINDOW:4:20 \
 MINLEN:50
echo Finish reads trimming
####### reads alignment #######
echo Begin reads alignment...
bowtie2-build --quiet --threads 8 $ref ref
bowtie2 --threads 8 --phred33 --no-discordant --no-unal --omit-sec-seq -x ref -1 FP.fq.gz -2 RP.fq.gz -S reads.sam
echo Finish reads alignment 
####### filter aligned reads #######
echo Begin aligned reads filtering...
samtools view -bST $ref -q 1 -F 0x0100 -@ 8 reads.sam > reads.bam
samtools sort -o reads.sort.bam -T samtmp -O bam -@ 8 reads.bam
rm reads.bam
bamtools filter -in reads.sort.bam -out reads.sort.filter1.bam -tag "XM:<3"
bamtools filter -in reads.sort.filter1.bam -out reads.sort.filter2.bam -tag "NM:<4" 
#XM: The number of mismatches in the alignment. 
#NM: The edit distance; that is, the minimal number of one-nucleotide edits (substitutions, insertions and deletions) needed to transform the read string into the reference string.
samtools stats -c 1,1000,1 reads.sort.filter2.bam > stats.txt
abn=`grep '^SN' stats.txt |cut -f2-|grep '^bases mapped (cigar)'|cut -f2`
avg=`grep '^COV' stats.txt |cut -f3-|sort -nrk2,2|head -1|cut -f1`
ulim=`expr $avg \\* 3`
blim=`expr $avg \\/ 3`
echo Finish aligned reads filtering
####### call SNPs #######
echo Begin SNP calling and filtering...
samtools mpileup -u -q 20 -Q 20 -f $ref reads.sort.filter2.bam | bcftools call --threads 8 -cv -p 0.05 > snp0.vcf
covflt="DP<$blim||DP>$ulim"
bcftools filter -e $covflt snp0.vcf  > snp1.vcf
echo Filter SNPs by $covflt...
bcftools filter -e 'DP4[0]>0 || DP4[1]>0' snp1.vcf > snp2.vcf
bcftools filter -i ' (GT="1|1" || GT="1/1") && AF1=1 ' snp2.vcf > snp3.vcf
echo Finish SNP calling and filtering
####### divergence filter and count SNPs ########
cd ..
echo Begin divergence filtering and SNPs counting...
perl QSBSA.pl -f=./tmp/snp3.vcf -v=off
echo Finish divergence filtering and SNPs counting

