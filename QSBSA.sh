#!/bin/bash
#SBATCH --job-name=QSBSA
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --mail-user=qiaoshan.lin@uconn.edu
#SBATCH -o QSBSA_%j.out
#SBATCH -e QSBSA_%j.err

module load Trimmomatic
module load bowtie2
module load samtools
module load bamtools
module load bcftools
module load bedtools
module load R/3.4.3
module load blast
module load minimap2
module load emboss

##### Functions

usage()
{
    echo "Program: qsbsa (Qiao-Shan Bulk Segregation Analysis)"
    echo "Version: 2.0"
    echo "Usage:"
    echo "    sh /path/to/script/QSBSA.sh [options] -r <ref.fa> [-1 <reads_1.fq> -2 <reads_2.fq> | -f <reads.txt>]"
    echo "Options:"
    echo "    -b | --bsa    To turn on the bulk segregation analysis with another genome to map to. Example: -b /home/CAM/qlin/SL9_Genome/version/v2.0/SL9_consensus.fa"
    echo "    -l | --library    To set the library for mutant SNPs comparing. Default: /home/CAM/qlin/BSA/LF10_vcf.txt"
    echo "    -s | --step	To restart with a certain step(trim/align/call-SNPs/filter/annotate). The results from previous run will be overwritten. Example: -s align"
    echo "    -g | --gtf 	To provide the gtf file for annotation of SNPs or set this option to be unavailable(n). Default: /home/CAM/qlin/resource/LF10/LF10g_v2.0.gtf"
    echo "    -c | --codingseq  To provide the coding sequence for annotation of SNPs in case no gtf file exists. Default: /home/CAM/qlin/resource/LF10/LF10T_v1.2.fa"
    echo "    -h | --help	To see the usage."
    echo "Notes:"
    echo "    If there are more than 2 reads files, save all paired reads files(including paths) to one txt file (one pair per line, seperated by a space) and use option -f."
}

trim()
{
for n in "${!array[@]}"
do
	java -jar $Trimmomatic PE \
	-phred33 \
	-threads 8 \
	${array[$n]} \
	$n\_FP.fq.gz $n\_FU.fq.gz \
	$n\_RP.fq.gz $n\_RU.fq.gz \
	LEADING:20 \
	TRAILING:20 \
	SLIDINGWINDOW:4:20 \
	MINLEN:50
done
echo Finish reads trimming
}

align()
{
bamarray=()
for n in "${!array[@]}"
do
	bowtie2-build --quiet --threads 8 $ref ref
	bowtie2 --threads 8 --phred33 --no-discordant --no-unal -x ref -1 $n\_FP.fq.gz -2 $n\_RP.fq.gz -U $n\_FU.fq.gz,$n\_RU.fq.gz -S $n\_reads.sam
	samtools view -bST $ref -@ 8 $n\_reads.sam > $n\_reads.bam
	bamarray=("${bamarray[@]}" $n"_reads.bam")
done
samtools merge reads.bam ${bamarray[@]}
samtools sort -@ 8 -o reads.sort.bam -T samtmp -O bam reads.bam
rm *_reads.sam *_reads.bam reads.bam
bamtools filter -in reads.sort.bam -out reads.sort.filter1.bam -tag "XM:<3"
samtools depth -aa -q 20 -Q 20 reads.sort.filter1.bam > depth.txt
echo Finish reads alignment and filtering
}

call()
{
samtools stats -c 1,1000,1 reads.sort.filter1.bam > stats.txt
abn=`grep '^SN' stats.txt |cut -f2-|grep '^bases mapped (cigar)'|cut -f2`
avg=`grep '^COV' stats.txt |cut -f3-|sort -nrk2,2|head -1|cut -f1`
ulim=`echo "$(($avg*25/10))"`
blim=`echo "$(($avg*10/25))"`
bcftools mpileup -Ou --threads 8 -Q 20 -q 20 -f $ref reads.sort.filter1.bam | bcftools call -Ou --threads 8 -mv | bcftools filter -e '%QUAL<20' > snp0.vcf
snpcount0=`grep -c -v '^#' snp0.vcf`
echo Round 0: $snpcount0 SNPs
#covflt="DP<$blim||DP>$ulim"
covflt="(DP4[0]+DP4[1]+DP4[2]+DP4[3])<$blim||(DP4[0]+DP4[1]+DP4[2]+DP4[3])>$ulim"
bcftools filter -e $covflt snp0.vcf  > snp1.vcf
echo Filter SNPs by coverage: $covflt
snpcount1=`grep -c -v '^#' snp1.vcf`
echo Round 1: $snpcount1 SNPs

frqflt="(DP4[2]+DP4[3])/(DP[0]+DP[1]+DP4[2]+DP4[3])<0.25"
bcftools filter -e $frqflt snp1.vcf > snp2.vcf
bcftools norm -f $ref snp2.vcf -Ov -o snp2.norm.vcf
mv snp2.norm.vcf snp2.vcf
echo Filter SNPs by alt frequency
snpcount2=`grep -c -v '^#' snp2.vcf`
echo Round 2: $snpcount2 SNPs

awk '{split($10, a, ":");split(a[2], b, ",");if (b[1]>40&&b[2]>20&&b[3]==0) print}' snp2.vcf > snp3.vcf
grep '^#' snp2.vcf |cat - snp3.vcf > tmp
mv tmp snp3.vcf
echo Filter SNPs by homozygousity
snpcount3=`grep -c -v '^#' snp3.vcf`
echo Round 3: $snpcount3 SNPs

echo Finish SNP calling and filtering
}

filter()
{
perl ~/scripts/QSBSA.pl -f=snp3.vcf -v=$lib -o=SNP
echo Finish SNPs filtering
}

annotate()
{
if [ "$gtf" != "n" ]; then
	awk '{print "$1~/"$1"/&&$4<="$2"&&$5>="$2"&&$3~/CDS|intron/{print "NR"\"\\\\t\"$0}"}' SNP.vcf > cmd.txt
	pos=()
	while read line; do   pos+=("$line");   done < cmd.txt
	snp=()
	while read line; do   snp+=("$line");   done < SNP.vcf
	for p in "${!pos[@]}"
	do
		awk "${pos[$p]}" $gtf |cut -f4-6,10 |tr "\n" "\t"> snptmp
		awk "${pos[$p]}" $gtf > snptmp2
#		feature=`awk "${pos[$p]}" $gtf |cut -f4`
		if [[ -s "snptmp" ]];then 
			snptmp=()
			while read line; do   snptmp+=("$line");   done < snptmp2
			for s in "${snptmp[@]}"
			do
				echo "$s"|cut -f4-6,8,10 |tr "\n" "\t" > tmp
				feature=`echo $s|awk '{print $4}'`
				strand=`echo $s|awk '{print $8}'`
				if [[ $feature =~ "CDS" && $strand =~ "+" ]];then
					startpos=`echo $s|awk '{print $5}'`
		        	        startframe=`echo $s|awk '{print $9}'`
        			        snppos=`echo "${snp[$p]}"|awk '{print $2}'`
        			        snpframe="$(( (($snppos-$startpos)%3-$startframe+3)%3 ))"
                			snpchr=`echo "${snp[$p]}"|awk '{print $1}'`
                			normstart=`echo "$(($snppos-$snpframe))"`
                			normend=`echo "$(($normstart+2))"`
               				samtools faidx $ref $snpchr:$normstart-$normend | tail -n +2 |awk -v frame="$snpframe" '{split($1,a,"");b=tolower(a[frame+1]);a[frame+1]=b;print a[1]a[2]a[3]}'|tr "\n" "\t"> codontmp
					echo "${snp[$p]}"|tr "\n" "\t"|cat - tmp codontmp |sed 's/$/\n/' >> SNP.annotate.vcf.tmp
				elif [[ $feature =~ "CDS" && $strand =~ "-" ]];then
					startpos=`echo $s|awk '{print $6}'`
                                        startframe=`echo $s|awk '{print $9}'`
                                        snppos=`echo "${snp[$p]}"|awk '{print $2}'`
                                        snpframe="$(( (($startpos-$snppos)%3-$startframe+3)%3 ))"
                                        snpchr=`echo "${snp[$p]}"|awk '{print $1}'`
                                        normstart=`echo "$(($snppos+$snpframe))"`
                                        normend=`echo "$(($normstart-2))"`
                                        samtools faidx $ref $snpchr:$normend-$normstart | tail -n +2 |perl -lpe 'tr/ATGCatgc/TACGTACG/;$_=reverse $_'|awk -v frame="$snpframe" '{split($1,a,"");b=tolower(a[frame+1]);a[frame+1]=b;print a[1]a[2]a[3]}'|tr "\n" "\t"> codontmp
                                        echo "${snp[$p]}"|tr "\n" "\t"|cat - tmp codontmp |sed 's/$/\n/' >> SNP.annotate.vcf.tmp
				else
					echo "${snp[$p]}"|tr "\n" "\t"|cat - tmp |sed 's/$/\n/' >> SNP.annotate.vcf.tmp
				fi
			done
		else
			echo "${snp[$p]}"|tr "\n" "\t"|cat - snptmp |sed 's/$/\n/' >> SNP.annotate.vcf.tmp
		fi
	done
	rm cmd.txt snptmp tmp snptmp2 codontmp
	mv SNP.annotate.vcf.tmp SNP.annotate.vcf
	cp SNP.annotate.vcf ../candidates.annotate.vcf
else
	awk '{if($2>1000)print $1":"$2-1000"-"$2+1000;else print $1":1-"$2+1000}' SNP.vcf > pos.txt
	pos=()
	while read line; do   pos+=("$line");   done < pos.txt
	for p in "${pos[@]}"
	do
	        samtools faidx $ref $p >> flanking.fa
	done
	rm pos.txt
	#select alignments that have 100% identity, SNP within alignment region; and print the SNP position at transcript to the last column of file
	/isg/shared/apps/blast/ncbi-blast-2.7.1+/bin/blastn -subject $cod -query flanking.fa -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send pident mismatch length sstrand" | awk '$9>=98{if($2==2001&&$3<=1001&&$4>=1001&&$12~/plus/)print $0"\t"1001-$3+$7;else if($2==2001&&$3<=1001&&$4>=1001&&$12~/minus/)print $0"\t"$7-1001+$3;else if($1~/:1-/&&$3<=($2-1000)&&$4>=($2-1000)&&$12~/plus/)print $0"\t"$2-1000-$3+$7;else if($1~/:1-/&&$3<=($2-1000)&&$4>=($2-1000)&&$12~/minus/)print $0"\t"$7-$2+1000+$3;else if($3<=1001&&$4>=1001&&$12~/plus/)print $0"\t"1001-$3+$7;else if($3<=1001&&$4>=1001&&$12~/minus/)print $0"\t"$7-1001+$3;}' > blast.txt
	cut -f5 blast.txt |sort|uniq > transcripts.txt
	perl ~/scripts/fetchSeq4fa2.pl $cod transcripts.txt blast.fa
	rm transcripts.txt flanking.fa
	cp blast.fa ../candidates.fa
fi
echo Finish annotation
}

bsa_annotate()
{
awk '$3>0{print $1"\t"$2"\t"$2+20000}' SNP.txt |tail -n +2 > tmp.bed
bedtools merge -i tmp.bed |awk '$3-$2>80000{print}'|bedtools merge -i - -d 200000 |awk '{print $1"\t"$2"\t"$3}'> high_count_regions.bed
rm tmp.bed
cd ..
awk '{print $1"\t"$2"\t"$2}' SNP.vcf |bedtools intersect -wb -a bsatmp/high_count_regions.bed -b - |awk '{print $4"\t"$5}' > candidates.txt
if [ "$gtf" != "n" ]; then
	array=()
	while read line; do   array+=("$line");   done < candidates.txt
	cp SNP.annotate.vcf candidates.annotate.vcf
	for n in "${array[@]}"
	do
		chr=`echo $n|tr " " "\t"|cut -f1`
		pos=`echo $n|tr " " "\t"|cut -f2`
		awk -v chr="$chr" -v pos="$pos" '{if($1 ~ chr&&$2 ~ pos){print $0"\thigh_count_region"}else{print}}' candidates.annotate.vcf > tmp.annotate.vcf
		mv tmp.annotate.vcf candidates.annotate.vcf 
	done
	mv candidates.annotate.vcf ../candidates.annotate.vcf
else
	awk '{if($2>1000)print $1":"$2-1000"-"$2+1000;else print $1":1-"$2+1000}' candidates.txt > pos.txt
        pos=()
        while read line; do   pos+=("$line");   done < pos.txt
        for p in "${pos[@]}"
    	do
            	samtools faidx $ref $p >> flanking.fa
    	done
    	rm pos.txt
    #select alignments that have 100% identity, SNP within alignment region; and print the SNP position at transcript to the last column of file
    	/isg/shared/apps/blast/ncbi-blast-2.7.1+/bin/blastn -subject $cod -query flanking.fa -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send pident mismatch length sstrand" | awk '$9>=98{if($2==2001&&$3<=1001&&$4>=1001&&$12~/plus/)print $0"\t"1001-$3+$7;else if($2==2001&&$3<=1001&&$4>=1001&&$12~/minus/)print $0"\t"$7-1001+$3;else if($1~/:1-/&&$3<=($2-1000)&&$4>=($2-1000)&&$12~/plus/)print $0"\t"$2-1000-$3+$7;else if($1~/:1-/&&$3<=($2-1000)&&$4>=($2-1000)&&$12~/minus/)print $0"\t"$7-$2+1000+$3;else if($3<=1001&&$4>=1001&&$12~/plus/)print $0"\t"1001-$3+$7;else if($3<=1001&&$4>=1001&&$12~/minus/)print $0"\t"$7-1001+$3;}' > blast.txt
    	cut -f5 blast.txt |sort|uniq > transcripts.txt
    	perl ~/scripts/fetchSeq4fa2.pl $cod transcripts.txt blast.fa
    	rm transcripts.txt flanking.fa
	cp blast.fa ../candidates.fa
fi
echo Finish annotation
}

##### Main

lib=/home/CAM/qlin/BSA/LF10_vcf.txt
gtf=/home/CAM/qlin/resource/LF10/LF10g_v2.0.gtf
cod=/home/CAM/qlin/resource/LF10/LF10T_v1.2.fa

if [ $# -lt 4 ]; then 
	usage
	exit 1
else
	while [ "$1" != "" ]; do
	    case $1 in
	        -r | --reference )      shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then 
						echo "missing value for -r"
						exit 1
					fi 
					ref=$1 
					;;
			-1 )                    shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then 
						echo "missing value for -1"
						exit 1
					fi
					read1=$1 
					;;
			-2 )                    shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then 
						echo "missing value for -2"
						exit 1
					fi
					read2=$1 
					;;
			-f | --file)		shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then 
						echo "missing value for -f"
						exit 1
					fi
					file=$1 
					;;
			-b | --bsa ) 	shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then
                        echo "missing value for -b"
                        exit 1
                    fi
                    bsaref=$1
                    ;;		
			-l | --library ) 	shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then
                        echo "missing value for -l"
                        exit 1
                    fi
                    lib=$1
                    ;;
 	       -s | --step )		shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then 
						echo "missing value for -s"
						exit 1
					fi
					step=$1 
					;;
			-g | --gtf ) 		shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then
						echo "missing value for -g"
						exit 1
					fi
					gtf=$1
					;; 
			-c | --codingseq ) 	shift
					if [[ "$1" == "" || "$1" =~ ^- ]]; then
						echo "missing value for -c"
						exit 1
					fi
					cod=$1
					;;
		    -h | --help )           usage
						exit 
					;;
	    	* )                     usage
						exit 1
	    esac
	    shift
	done
fi

#if [[ ( ( "$file" == "" && "$read1" == "") || ("$file" != "" && "$read1" != "") ) || ( ( "$read1" != "" && "$read2" == "" ) || ( "$read2" != "" && "$read1" == "" ) ) ]]; then
#	usage
#	exit 1
#fi

if [ "$file" != "" ]; then
	sed -i '/^$/d' $file
	array=()
	while read line; do   array+=("$line");   done < $file
else 
	array=( "$read1 $read2" )
fi 

if [[ `ls -l | grep tmp` ]]; then echo overwrite previous results;else mkdir tmp;fi
cd tmp

case $step in
        "" ) ;&
        "trim" ) trim ;&
        "align" ) align ;&
        "call-SNPs" ) 
		call 
		newvcf=`pwd`
		echo "$newvcf/snp2.vcf">> /home/CAM/qlin/BSA/vcf.txt
		sort /home/CAM/qlin/BSA/vcf.txt |uniq > /home/CAM/qlin/BSA/newvcf.txt
		mv /home/CAM/qlin/BSA/newvcf.txt /home/CAM/qlin/BSA/vcf.txt
		echo VCF library updated at /home/CAM/qlin/BSA/vcf.txt
		;&
        "filter" ) filter ;&
	"annotate" ) annotate ;;
        *) echo "Error: no such a step"; exit 1 ;;
esac

if [ "$bsaref" != "" ]; then
	mkdir bsatmp
	cd bsatmp
	echo Begin bulk segregation analysis......
	oldref=$ref
	ref=$bsaref
	case $step in
        "" ) ;&
        "trim" ) cp ../*fq.gz ./ ;&
        "align" ) align ;&
        "call-SNPs" ) call ;&
        "filter" ) 
		perl ~/scripts/QSBSA.pl -f=snp3.vcf -m=on -v=off -o=SNP 
		cp SNP.eps ../../SNP.eps
		echo Finish SNPs filtering
		;&
	"annotate" ) bsa_annotate ;;
        *) echo "Error: no such a step"; exit 1 ;;
	esac
fi

grep 'intron' ../candidates.annotate.vcf |awk '$2-$12<5||$13-$2<5{print}' |cut -f1-2,4-5,11- > ../check_first.vcf
grep 'CDS' ../candidates.annotate.vcf |cut -f1-2,4-5,11- >> ../check_first.vcf


