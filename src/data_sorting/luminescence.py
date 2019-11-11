#!/usr/bin/env python
# coding: utf-8

# <b> Sort the protoplast luminescence data from the xlsx output from the Glariostar platereader. 

# In[ ]:


import pandas as pd
import os
from pathlib import Path
import glob


# In[1]:





# In[83]:


def xlsx_2_csv(xlsx):  
    """ Function to read and convert xlsx file to csv file """
    
    # Read in the xlsx file, second sheet
    file = pd.read_excel(xlsx, 'End point', index_col=None) 
    
    filename = os.path.basename(xlsx)
    removed_extension = os.path.splitext(filename)[0]
    parent = Path(xlsx).parent.parent #find parent directory to the one the xlsx fiels are in
    
    file.to_csv(f'{parent}/csvs/{removed_extension}.csv', encoding='utf-8', index=False)


# In[84]:


#find all xlsx files recursively in the 'to_be_sorted' folder
xlsx_filenames = glob.glob('/home/witham/Documents/Pipeline/data/luminescence/to_be_sorted/**/*.xlsx', recursive=True)


# In[85]:


#run the xlsx_2_csv function across all xlsx file in to_be_sorted folder
list(map(xlsx_2_csv,xlsx_filenames))               
                 


# use os.scandir when scanning a directory, this is the fastest way according to Matt

# In[ ]:




