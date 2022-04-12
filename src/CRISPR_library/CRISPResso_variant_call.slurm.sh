#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -t 1-08:0
#SBATCH -c 32
#SBATCH --mem=30G
#SBATCH -J pacbio_variantcall
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=samuel.witham@earlham.ac.uk
#SBATCH -o pacbio_variantcall
#SBATCH -e  /hpc-home/witham/pacbio/data/CRISPR_library/pacbio/demultiplexed/Data_Package_Batch_04_04_2022/log/pacbio_variantcall.err


source /ei/projects/1/1ee4c3b1-03f3-4307-a995-693a8bcf9888/software/miniconda3/bin/activate

bash /hpc-home/witham/pacbio/src/CRISPR_library/./pacbio_variantcall_hpc_test.sh
