#!/bin/bash
#SBATCH -p general
#SBATCH -J homer
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 0-10:00
#SBATCH --mem 32000
#SBATCH -o homer.out
#SBATCH -e homer.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dsong@hsph.harvard.edu

source new-modules.sh
module load homer/latest-fasrc01

#sed '/GL/d' ~/BST281/homer/H3K27ac_increase.bed > ./newfile.bed
#awk '$0="chr"$0' newfile.bed > ./newfile1.bed
#awk 'BEGIN {OFS="\t"}{$6=0; print $0}' newfile1.bed >~/BST281/homer/H3K27ac_increase_cl.bed
#rm -f newfile.bed newfile1.bed
#findMotifsGenome.pl ~/BST281/homer/H3K27ac_increase_cl.bed hg19 ~/BST281/homer_output/H3K27ac_increase -preparsedDir ~/BST281/homer/preparsed -size 500

#sed '/GL/d' ~/BST281/homer/H3K4me3_increase.bed > ./newfile.bed
#awk '$0="chr"$0' newfile.bed > ./newfile1.bed
#awk 'BEGIN {OFS="\t"}{$6=0; print $0}' newfile1.bed >~/BST281/homer/H3K4me3_increase_cl.bed
#rm -f newfile.bed newfile1.bed
#findMotifsGenome.pl ~/BST281/homer/H3K4me3_increase_cl.bed hg19 ~/BST281/homer_output/H3K4me3_increase -preparsedDir ~/BST281/homer/preparsed -size 500

#sed '/GL/d' ~/BST281/homer/H3K27ac_decrease.bed > ./newfile.bed
#awk '$0="chr"$0' newfile.bed > ./newfile1.bed
#awk 'BEGIN {OFS="\t"}{$6=0; print $0}' newfile1.bed >~/BST281/homer/H3K27ac_decrease_cl.bed
rm -f newfile.bed newfile1.bed
findMotifsGenome.pl ~/BST281/homer/H3K27ac_decrease_cl.bed hg19 ~/BST281/homer_output/H3K27ac_decrease -preparsedDir ~/BST281/homer/preparsed -size 500

sed '/GL/d' ~/BST281/homer/H3K4me3_decrease.bed > ./newfile.bed
awk '$0="chr"$0' newfile.bed > ./newfile1.bed
awk 'BEGIN {OFS="\t"}{$6=0; print $0}' newfile1.bed >~/BST281/homer/H3K4me3_decrease_cl.bed
rm -f newfile.bed newfile1.bed
findMotifsGenome.pl ~/BST281/homer/H3K4me3_decrease_cl.bed hg19 ~/BST281/homer_output/H3K4me3_decrease -preparsedDir ~/BST281/homer/preparsed -size 500
