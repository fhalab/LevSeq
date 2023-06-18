import os 
import numpy as np
#import holoviews as hv
#hv.extension('bokeh')
#import pandas as pd
import re
import subprocess
import glob

# Activate medaka
#subprocess.call('conda run -n medaka', shell = True)

# Obtain path to current basecalled fasta file folder
crr_pwd = os.getcwd()
subprocess.call('mkdir consensus',shell = True)
# Check if fast5 folder exist
if os.path.exists(crr_pwd + '/fast5') == True:
    os.chdir('fast5')
    crr_pwd = os.getcwd()
else:
    crr_pwd = crr_pwd

# Obtain path to fastq files that passed basecalling
os.chdir('pass')
pass_pwd = os.getcwd()

# Obtain path to front barcode folder
os.chdir('fbc')
fbc_pwd = os.getcwd()

# Read processed forward and reverse barcode folder
candidates_all = glob.glob(os.path.join('barcode*/barcode*'))
# Loop through each folder and create consensus sequences
for candidate in candidates_all:
    fbc = candidate.split(os.sep)[0][-2::]
    rbc = candidate.split(os.sep)[1][-2::]
    bc_pair = os.path.join(candidate, fbc + rbc  + ".fastq")

# Set parameters for medaka consensus
    mc_input = bc_pair
    mc_output = fbc_pwd + f'/{candidate}'
    mc_temp = '~/tam-lq.fasta'
    mc_model = 'r1041_e82_400bps_hac_g615'

# Call medaka consensus
    subprocess.call(f'medaka_consensus -i {mc_input} -d {mc_temp} -o {mc_output} -m {mc_model} -r -f -x',shell = True)


# Set parameters for medaka stitch
    ms_input = fbc_pwd + '/' +  candidate + '/consensus_probs.hdf'
    ms_output = crr_pwd + '/consensus/' + fbc + rbc + '.fasta'
    ms_temp = '~/tam-lq.fasta'

# Call medaka stitch
    subprocess.call(f'medaka stitch {ms_input} {ms_temp} {ms_output} --qualities', shell = True)
