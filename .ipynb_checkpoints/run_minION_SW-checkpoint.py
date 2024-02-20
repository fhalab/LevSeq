#!/usr/bin/env python
# coding: utf-8

# # Every Variant Sequencing with Oxford Nanopore Technologies
# 
# This script is being used after sequencing. The raw pod5 files can be basecalled or the already basecalled files can be used directly (fastq.gz)
# 
# ## Workflow
# 
# ### 1. Basecalling (Optional)
# 
# - The raw reads are stored in the main folder of ONT (e.g /var/lib/minknow/data). Enter the experiment name as input. 
# - Sequences are basecalled based on the model of choice. If enough computational power is available, we recommend "sup" method
# 
# ### 2. Demultiplexing 
# - Each reead is assigned to a well/plate combination. 
# 
# ### 3. Variant Calling
# - Minimap2 for creating Multiple Sequence Alignment (MSA)
# - Base Frequency Caller is being used for variant calling
# 
# 

# ### Packages 

# In[45]:


# Import all packages

import sys
sys.path.append("/home/longy/git/MinION")

from minION.util import IO_processor
from minION.basecaller import Basecaller

from minION.variantcaller import *

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import importlib
importlib.reload(IO_processor)


# ### Meta Data 
# 
# - Provide the following arguments:
# 
# - Result Path: Path where the minion result folder will be created. All experiment results are then stored within the folder
# - Experiment Name: The experiment name is assigned when running the sequencer. Use the same name for identification
# 

# In[46]:


result_path = Path("/home/longy/")
experiment_name = "20240126-RL-8-63"
basecall_model_type = "sup"
result_folder = IO_processor.create_folder( experiment_name,
                                            basecall_model_type, 
                                            target_path=result_path)




# Create Barcode fasta file 
barcode_path = "minION/barcoding/minion_barcodes.fasta" # Path to standard barcode file
front_prefix = "NB"
back_prefix = "RB"
bp = IO_processor.BarcodeProcessor(barcode_path, front_prefix, back_prefix)
#barcode_path = result_folder / "minion_barcodes_filtered.fasta"

# Barcode indexes
front_min = 1
front_max = 96
back_min = 9
back_max = 12

bp.filter_barcodes(barcode_path, (front_min,front_max), (back_min,back_max))


file_to_experiment= f"/var/lib/minknow/data/{experiment_name}"
template_fasta = "/home/emre/PgA9.fasta"

# Basecalling
basecall_folder = result_folder / "basecalled"
basecall_folder.mkdir(parents=True, exist_ok=True)
experiment_folder = IO_processor.find_experiment_folder(experiment_name) # Folder where pod5 files are located

# Demultiplexing
experiment_name = experiment_name + "_" + basecall_model_type
result_folder_path = IO_processor.find_folder(result_path, experiment_name)


# In[47]:


# Add conditions to avoid running the script accidentally
skip_basecalling = True
skip_demultiplex = True
skip_variant_calling = False


# In[48]:


result_folder


# ### Step 1 (Optional): Basecall reads
# 
# - Basecall can usually be done while sequencing (if GPU available?)
# - Otherwise, basecall afterwards

# In[50]:


if not skip_basecalling:
    pod5_files = IO_processor.find_folder(experiment_folder, "pod5")
    bc = Basecaller(basecall_model_type, pod5_files, basecall_folder, fastq = True)
    bc.run_basecaller()


# In[51]:


# Find fastq files
file_to_fastq = IO_processor.find_folder(experiment_folder, "fastq_pass")
print(file_to_fastq)


# ### Step 2: Demultiplex with SW
# - Demultiplex with SW 

# In[52]:


if not skip_demultiplex:
    path_to_code = "/home/emre/github_repo/MinION/source/source/demultiplex"
    prompt = f"{path_to_code} -f {file_to_fastq} -d {result_folder} -b {barcode_path} -w {100} -r {100}"
    subprocess.run(prompt, shell=True)


# In[53]:


demultiplex_folder = result_folder 
print(demultiplex_folder)


# ### Step 3: Call Variant with PileUP Analysis
# 
# - Call Variant with min freq of 0.4 & depth min 15

# Read Summary file (Optional):
# 

# In[54]:


demultiplex_folder_name = result_folder


# In[55]:


if not skip_variant_calling:
    vc = VariantCaller(experiment_folder, 
                   template_fasta, 
                   demultiplex_folder_name=demultiplex_folder_name, 
                   padding_start=0, 
                   padding_end=0)
    
    variant_df = vc.get_variant_df(qualities=True, 
                                threshold=0.2,
                                min_depth=5)
    seq_gen = IO_processor.SequenceGenerator(variant_df, template_fasta)
    variant_df = seq_gen.get_sequences()
    #TODO: Save the variant_df to a file after running. Currently it is not saved.


# In[56]:


#20 - 30
variant_df.to_csv(result_folder / "variant_df.csv", index=False)  


# In[57]:


variant_df


# In[ ]:




