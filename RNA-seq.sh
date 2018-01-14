#!/bin/bash
input_path=$1
fastq_dump=/gluster/home/yuanxiao/software/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump.2.8.0
fastqc=/gluster/Bioshare/Biotools/QC/FastQC/fastqc
qc_stat=/gluster/Bioshare/Biotools/QC/fastq_stat
adapter=/gluster/home/yuanxiao/software/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa
trim=/gluster/home/yuanxiao/software/Trimmomatic-0.35/trimmomatic-0.35.jar

if [ "$2" = "hg38" ];then
gtf=/gluster/home/yuanxiao/other/Homo_sapiens.GRCh38.87.chr+chr.gtf
index=/gluster/home/yuanxiao/ref/hg38_for_HISAT/genome
cd ${input_path}
echo -e "gtf="$gtf "\nindex="$index >ref.log

elif [ "$2" = "hg19" ];then
gtf=/gluster/home/yuanxiao/other/Homo_sapiens.GRCh37.75+chr.gtf
index=/gluster/home/yuanxiao/ref/hg19_for_HISAT/genome
cd ${input_path}
echo -e "gtf="$gtf "\nindex="$index >ref.log

elif [ "$2" = "mm10" ];then
gtf=/gluster/home/yuanxiao/other/Mus_musculus.GRCm38.83+chr.gtf
index=/gluster/home/yuanxiao/ref/mm10_for_HISAT/genome
cd ${input_path}
echo -e "gtf="$gtf "\nindex="$index >ref.log

else
cd ${input_path}
echo "error" > error.log
fi





#Step 0:Decompress
: << !
for i in $(ls ${input_path}/0.data/*.sra)
do
$fastq_dump --split-3 $i -O ${input_path}/0.data
done
!

for i in $(ls ${input_path}/0.data/*.gz)
do
sample_name=`basename $i|awk -F"." '{print $1}'`
gzip -cd $i >>${input_path}/0.data/${sample_name}.fastq
done






#Step 1: Trim
mkdir ${input_path}/1.Trim
mkdir ${input_path}/1.Trim/before_trim
mkdir ${input_path}/1.Trim/after_trim
for i in $(ls ${input_path}/0.data/*.fastq)
do
$fastqc $i -t 10 -o ${input_path}/1.Trim/before_trim
done

if [ "$3" = "SE" ];then
for i in $(ls ${input_path}/0.data/*.fastq)
do
sample_name=`basename $i|sed 's/.fastq//g'`
java -jar $trim SE -threads 24 ${input_path}/0.data/${sample_name}.fastq ${input_path}/1.Trim/${sample_name}_trimmed ILLUMINACLIP:$adapter:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
$fastqc ${input_path}/1.Trim/${sample_name}_trimmed -t 10 -o ${input_path}/1.Trim/after_trim
#$qc_stat 
done

elif [ "$3" = "PE" ];then
for i in $(ls ${input_path}/0.data/*_R1.fastq)
do
sample_name=`basename $i|sed 's/_R1.fastq//g'`
java -jar $trim PE -threads 24 ${input_path}/0.data/${sample_name}_R1.fastq ${input_path}/0.data/${sample_name}_R2.fastq -baseout ${input_path}/1.Trim/${sample_name}_trimmed ILLUMINACLIP:$adapter:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
$fastqc ${input_path}/1.Trim/${sample_name}_trimmed_1P -t 10 -o ${input_path}/1.Trim/after_trim
$fastqc ${input_path}/1.Trim/${sample_name}_trimmed_2P -t 10 -o ${input_path}/1.Trim/after_trim
$qc_stat ${input_path}/1.Trim/${sample_name}_trimmed_1P ${input_path}/1.Trim/${sample_name}_trimmed_2P ${input_path}/1.Trim/${sample_name}_trimmed_stat
done

else
cd ${input_path}
echo "error" > error.log
fi
Rscript /gluster/home/yuanxiao/scripts/R/QC_summary.R ${input_path}/1.Trim





#Step 2: HISAT2: Pay attention to strand-specific sequence!!! --rna-strandness FR RF
mkdir ${input_path}/2.HISAT2
if [ "$3" = "SE" ];then
for i in $(ls ${input_path}/1.Trim/*_trimmed)
do
sample_name=`basename $i|sed 's/_trimmed//'`
hisat2 -p 24 --dta -x $index -U $i -S ${input_path}/2.HISAT2/${sample_name}.sam 2>${input_path}/2.HISAT2/${sample_name}.uniq_mapping_rate
samtools sort -@ 24 -o ${input_path}/2.HISAT2/${sample_name}.bam ${input_path}/2.HISAT2/${sample_name}.sam
done

elif [ "$3" = "PE" ];then
for i in $(ls ${input_path}/1.Trim/*_trimmed_1P)
do
sample_name=`basename $i|sed 's/_trimmed_1P//g'`
hisat2 -p 24 --dta -x $index -1 ${input_path}/1.Trim/${sample_name}_trimmed_1P -2 ${input_path}/1.Trim/${sample_name}_trimmed_2P -S ${input_path}/2.HISAT2/${sample_name}.sam 2>${input_path}/2.HISAT2/${sample_name}.uniq_mapping_rate
samtools sort -@ 24 -o ${input_path}/2.HISAT2/${sample_name}.bam ${input_path}/2.HISAT2/${sample_name}.sam
done

else
cd ${input_path}
echo "error" > error.log
fi





#Step 3: StringTie
mkdir ${input_path}/3.StringTie
for i in $(ls ${input_path}/2.HISAT2/*.bam)
do
sample_name=`basename $i|sed 's/.bam//'`
samtools flagstat $i >${input_path}/2.HISAT2/${sample_name}.mapping_rate.stat
mkdir ${input_path}/3.StringTie/${sample_name}
stringtie -e -B -p 24 -G $gtf -o ${input_path}/3.StringTie/${sample_name}/${sample_name}.gtf $i
done

Rscript /gluster/home/yuanxiao/scripts/R/mapping_rate_summary.R ${input_path}/2.HISAT2
Rscript /gluster/home/yuanxiao/scripts/R/mapping_rate_summary2.R ${input_path}/2.HISAT2
python /gluster/home/yuanxiao/scripts/python/prepDE.py -i ${input_path}/3.StringTie -g ${input_path}/gene_count_matrix.csv -t ${input_path}/transcript_count_matrix.csv





#Step 4: DESeq2
#Rscript /gluster/home/yuanxiao/scripts/RNA-seq_pipeline_1/deseq2.R ${input_path}




