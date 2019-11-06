#!/bin/bash -e
#SBATCH -p ei-short
#SBATCH -t 0-00:03
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -J x_promz
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=samuel.witham@earlham.ac.uk
#SBATCH -o x_promz
#SBATCH -e  /ei/workarea/group-eg/project_PromoterArchitecturePipeline/log/x_promz.err


source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate
conda activate PromoterArchitecturePipeline
python /ei/workarea/group-eg/project_PromoterArchitecturePipeline/scripts/extract_promoter.py
