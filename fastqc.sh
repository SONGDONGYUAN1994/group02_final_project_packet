#!/bin/bash
#SBATCH -p general
#SBATCH -J fastqc
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 0-10:00
#SBATCH --mem 8000
#SBATCH -o fastqc.out
#SBATCH -e fastqc.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dsong@hsph.harvard.edu

cd /n/home08/songdongyuan/BST281/

source new-modules.sh
module load fastqc/0.11.5-fasrc01

fastqc -o ~/BST281/fastqc_output -t 16 ~/BST281/chip/Alignment_Post_Processing_15005.bam ~/BST281/chip/Alignment_Post_Processing_15009.bam ~/BST281/chip/Alignment_Post_Processing_15022.bam ~/BST281/chip/Alignment_Post_Processing_15175.bam ~/BST281/chip/Alignment_Post_Processing_15180.bam ~/BST281/chip/Alignment_Post_Processing_15193.bam ~/BST281/chip/Alignment_Post_Processing_15223.bam ~/BST281/chip/Alignment_Post_Processing_15280.bam

