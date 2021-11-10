#!/bin/bash


# alignment for pacbio reads 
# param:
# 1: reads file
# 2: prefix 
# 3: reference (hg38)
function align_pacbio()
{
    # ~/minimap/minimap2 -ax map-pb $3 $1 > ${2}.sam
    # samtools view -F 2308 -b -T $3 ${2}.sam > ${2}.bam
    ngmlr -t 8 -r $3 -q $1 -o ${2}_pb.sam
    samtools sort -O BAM ${2}_pb.sam -o ${2}_pb.bam --threads 8
    samtools index ${2}_pb.bam

}

# alignment for nanopore reads 
# param:
# 1: reads file
# 2: prefix 
# 3: reference (hg38)
function align_nanopore()
{
    # /home/eric/minimap/minimap2 -ax map-ont $3 $1 > ${2}.sam
    # samtools view -F 2308 -b -T $3 ${2}.sam > ${2}.bam
    # rm ${2}.sam
    ngmlr -t 8 -r $3 -q $1 -o ${2}_ont.sam -x ont
    samtools sort -O BAM ${2}_ont.sam -o ${2}_ont.bam --threads 8
    samtools index ${2}_ont.bam
}

# alignment for hic paired-end reads 
# param:
# 1: forward reads
# 2: reverse reads
# 3: prefix 
# 4: reference (hg19)
function align_hic()
{
    bwa mem -t 8 -5SP $4 $1 $2 | samtools view -F 2304 -b -T $4 | samtools sort --thread 8 -n | samtools fixmate -O bam - ${3}.fixmate.bam
    java -jar picard.jar MarkDuplicates I=./${3}.fixmate.bam O=./${3}.bam M=./${3}.markdup.txt
    rm  ./${3}.fixmate.bam
}

# alignment for ngs paired-end reads 
# param:
# 1: forward reads
# 2: reverse reads
# 3: prefix 
# 4: reference
function align_ngs()
{
    # bwa mem -t 8 -Y -M $4 $1 $2 | samtools view -F 2304 -b -T $4 | samtools sort --thread 8 -n | samtools fixmate -O bam - ${3}.fixmate.bam
    # java -jar picard.jar MarkDuplicates I=./${3}.fixmate.bam O=./${3}.bam M=./${3}.markdup.txt
    # rm  ./${3}.fixmate.bam
    bwa mem -t 8 $4 $1 $2 > ${3}.sam
    samtools sort -O BAM ${3}.sam -o ${3}.bam --threads 8
    samtools index ${3}.bam
}

# subsampling for bam file
# 1: subsampling fraction
# 2: input bam file
# 3: prefix
function subsampling()
{
    samtools view -s $1 -b $2 > ${3}.sam
    samtools sort -O BAM ${3}.sam -o ${3}.bam --threads 8
    samtools index ${3}.bam
}

# mix normal data with tumor data
# 1: normal BAM
# 2: normal purity
# 3: tumor BAM
# 4: tumor purity
# 5: prefix

function merge()
{
    samtools  view -b -s $2 $1 > temp1.bam
    samtools  view -b -s $4 $3 > temp2.bam
    samtools merge ${5}.bam temp1.bam temp2.bam -f
    samtools index ${5}.bam
    rm -f temp1.bam
    rm -f temp2.bam
}

# alignment for ngs paired-end reads 
# param:
# 1: fastq_dir
# 2: reference_dir
# 3: prefix 
# 4: gatk path
function align_tenx()
{
    longranger wgs \
    --id=$3 \
    --sample=$3 \
    --fastqs=$1 \
    --reference=$2 \
    --vcmode=$4 \
    --jobmode=local \
    --localcores=24 \
    --localmem=96
}
# PE
# align_ngs ./bfb_r1.fq.gz ./bfb_r2.fq.gz bfb ./hg38.fa

# coverage
# subsampling 0.67 ./bfb.bam bfb_20x
# subsampling 0.33 ./bfb.bam bfb_10x
# subsampling 0.17 ./bfb.bam bfb_5x

# tumor purity
# align_ngs ./normal_r1.fq.gz ./normal_r2.fq.gz normal ./hg38.fa
# merge ./normal.bam 0.25 ./bfb.bam 0.75 bfb_75x
# merge ./normal.bam 0.5 ./bfb.bam 0.5 bfb_50x
# merge ./normal.bam 0.75 ./bfb.bam 0.25 bfb_25x
# merge ./normal.bam 0.8 ./bfb.bam 0.2 bfb_20x
# merge ./normal.bam 0.9 ./bfb.bam 0.1 bfb_10x

# PB
# 1: reads file
# 2: prefix 
# 3: reference (hg38)
# align_pacbio ./bfb_pb_0001.fastq bfb ./hg38.fa

# ONT
# 1: reads file
# 2: prefix 
# 3: reference (hg38)
align_nanopore ./bfb1_ont/pass.fastq bfb1_ont ./hg38.fa