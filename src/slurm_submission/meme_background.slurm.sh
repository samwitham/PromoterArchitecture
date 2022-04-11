#!/bin/bash -e
#SBATCH -p ei-short
#SBATCH -t 0-00:03
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -J x_promz
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=sam.witham@earlham.ac.uk
#SBATCH -o x_promz
#SBATCH -e  /ei/workarea/group-eg/project_PromoterArchitecturePipeline/log/background.err


source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate
conda activate PromoterArchitecturePipeline
bash /ei/workarea/group-eg/project_PromoterArchitecturePipeline/scripts/meme_suite/./background.sh /ei/workarea/group-eg/project_PromoterArchitecturePipeline/data/genomes/promoters.gff3
