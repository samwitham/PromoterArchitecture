#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -t 0-08:0
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH -J whole_pipeline
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=samuel.witham@earlham.ac.uk
#SBATCH -o whole_pipeline
#SBATCH -e  /ei/workarea/group-eg/project_PromoterArchitecturePipeline/log/whole_pipeline.err


source /ei/workarea/group-eg/project_PromoterArchitecturePipeline/software/miniconda3/bin/activate

bash /ei/workarea/group-eg/project_PromoterArchitecturePipeline/src/whole_pipeline/./whole_pipeline.sh
