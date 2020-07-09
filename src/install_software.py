#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#pip install pybiomart


# In[1]:


#conda install -y -c future requests requests-cache

conda create --name PromoterArchitecturePipeline python=3
conda activate PromoterArchitecturePipeline
conda install -c anaconda mysql git biopython pymysql configargparse 
conda install -c bioconda gffutils 

conda install --channel bioconda bedtools htslib bedops pybedtools bcbiogff pyfaidx

conda install -c conda-forge git-lfs
#install statannot
pip install statannot
#in project_PromoterArchitecturePipeline folder:

# git lfs install
# git lfs track "data/genomes"
# git pull
# git status
# git config --global credential.helper 'cache --timeout=3600'
# git pull -u origin master
# git add .gitattributes



