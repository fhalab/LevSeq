import glob
import numpy as np
import pandas as pd
import re
import os
import subprocess
import holoviews as hv
hv.extension('bokeh')

# Obtain path to where the basecalled fasta file is stored
crr_pwd = os.getcwd()

# Basecall the reads
# Set parameters (potential user input)
input_file = crr_pwd
save_path = crr_pwd
config = 'dna_r10.4.1_e8.2_260bps_hac.cfg'
# Combine input variables into one string
cmd_basecall = f'guppy_basecaller -i {input_file} -s {save_path} -c {config} -x cuda:0,1:80% --chunks_per_runner 64'
# Call shell process guppy_basecaller
#subprocess.call(cmd_basecall, shell = True)

# Obtain path to passed basecalls
os.chdir('pass')
pass_pwd = os.getcwd()

# Identification of forward barcodes
# Set parameters for barcoding (potential user input)
barcode_kit = 'MY-CUSTOM-BARCODES'
fbc_input = pass_pwd
fbc_output = pass_pwd + '/fbc'
# Combine input forward barcode demultiplex variables into one string
cmd_fbc = f'guppy_barcoder --input_path {fbc_input} --save_path {fbc_output} --data_path ~/mydata/barcoding --barcode_kits {barcode_kit}'
# Call shell process guppy_barcoder
#subprocess.call(cmd_fbc, shell = True)

# Read processed forward barcode folders
os.chdir('fbc')
fbc_pwd = os.getcwd()
candidates_fw = glob.glob(os.path.join("barcode*"))
for candidate_fw in candidates_fw:
    fw = os.path.join(candidate_fw, candidate_fw[-2::] + ".fastq")
    # Concatenate all fastq files into 1
    with open(fw, 'w') as output_file: 
        fastq_fw = glob.glob(os.path.join(candidate_fw, "fastq_runid*"))
        for f in fastq_fw:
            with open(f, 'r') as input_file: 
                contents = input_file.read()	
                output_file.write(contents)

    # Identification of reverse barcodes
    # Set parameters for reverse barcodes(potential user input)
    rbc_input = fbc_pwd + f'/{candidate_fw}'
    rbc_output = rbc_input
    rbc_kit = 'MY-CUSTOM-BARCODES-REV'
    cmd_rbc = f'guppy_barcoder --input_path {rbc_input} --save_path {rbc_output} --data_path ~/mydata/barcoding --barcode_kits {rbc_kit}'
#    subprocess.call(cmd_rbc, shell = True)

# Create empty datafrane
df_count  = pd.DataFrame({'FBC':[],'RBC':[], 'Aligned':[]})
# Obtain forward and reverse barcode pair
candidates = glob.glob(os.path.join("barcode*/barcode*"))
for candidate in candidates:
    fbc = candidate.split(os.sep)[0][-2::]
    rbc = candidate.split(os.sep)[1][-2::]
    bc_pair = os.path.join(candidate, fbc + rbc  + ".fastq")

    # Concatenate all fastq files into 1
    with open(bc_pair,'w')as output_pair:
        fastq_pair = glob.glob(os.path.join(candidate, 'fastq_runid*'))
        for fp in fastq_pair:
            with open(fp, 'r') as input_pair:
                content_pair = input_pair.read()
                output_pair.write(content_pair)

    # Use regex to find and count sequences
    regex_ct = r'runid'
    f_pair = open(bc_pair, 'r')
    ct_seq = len(re.findall(regex_ct,f_pair.read()))

    # Append sequence count to dataframe
    data_temp = {'RBC':[rbc],'FBC':[fbc], 'Aligned':[ct_seq]}
    df_temp = pd.DataFrame(data_temp)
    df_count = pd.concat([df_temp,df_count], axis=0)

list_ct = df_count.values.tolist()
tp_ct = [tuple(x) for x in list_ct]
hm = hv.HeatMap(tp_ct).sort().opts(width = 600, invert_yaxis = True,xaxis = 'top')    
hm_ct = hm * hv.Labels(hm).opts(padding = 0)
hv.save(hm_ct, 'alignment_count.html')
