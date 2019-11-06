#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd


# In[4]:


def filter_fimo_qvalue(input_file,output_file,qvalue_cutoff):
    """function to read in FIMO .tsv file, filter values where q is less than the cutoffvalue and output a file"""
    df = pd.read_csv(input_file,sep='\t',header=0)
    print(df)
    
    
    
    


# In[5]:


filter_fimo_qvalue('/home/witham/Documents/ugene/Reference_sequences/N-responsive/fastas/fimo_out/fimo.tsv','/home/witham/Documents/ugene/Reference_sequences/N-responsive/fastas/fimo_out/fimo_filtered.tsv',0.05)


# In[ ]:




