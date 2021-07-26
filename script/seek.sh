#!/bin/bash
#SBATCH --partition=cpu_14d1n
#SBATCH --job-name="seeksv"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1
#SBATCH --output=seek%j.log
#SBATCH --qos=normal
samtools index -@ 48 $2
seeksv getclip -o $1 $2
bwa mem -t 128 $3  $1.clip.fq.gz | samtools view -@ 48 -Sb -o $1.clip.bam -
seeksv getsv $1.clip.bam $2 $1.clip.gz $1.seek.sv.txt $1.unmapped.clip.fq.gz
