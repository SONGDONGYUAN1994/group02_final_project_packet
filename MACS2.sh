#!/bin/bash
#SBATCH -p general
#SBATCH -J macs2
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -t 0-10:00
#SBATCH --mem 8000
#SBATCH -o macs2.out
#SBATCH -e macs2.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dsong@hsph.harvard.edu

cd /n/home08/songdongyuan/BST281/

source new-modules.sh
module load macs2/2.1.1.20160309-fasrc01

macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15005.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15005 -q 0.01
macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15009.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15009 -q 0.01
macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15022.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15022 -q 0.01
macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15175.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15175 -q 0.01
macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15180.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15180 -q 0.01
macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15193.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15193 -q 0.01
macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15223.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15223 -q 0.01
macs2 callpeak -t ~/BST281/chip/Alignment_Post_Processing_15280.bam --outdir ~/BST281/macs2_output -f BAMPE -g hs -n 15280 -q 0.01