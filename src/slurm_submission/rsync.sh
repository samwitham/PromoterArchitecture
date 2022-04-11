#!/bin/bash -e
#SBATCH -p ei-short
#SBATCH -t 0-00:10
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -J rsync
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=sam.witham@earlham.ac.uk
#SBATCH -o rsync
#SBATCH -e  /ei/projects/1/1ee4c3b1-03f3-4307-a995-693a8bcf9888/log/rsync.err


source /ei/projects/1/1ee4c3b1-03f3-4307-a995-693a8bcf9888/software/miniconda3/bin/activate
conda activate pacbio
rsync -a /ei/data/reads/PIP-2918/Data_Package_Batch_15_03_2022 /ei/projects/1/1ee4c3b1-03f3-4307-a995-693a8bcf9888/data/raw