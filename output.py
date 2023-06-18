import pandas as pd
import holoviews as hv
hv.extension('bokeh')
import glob
import re
import subprocess
import os
import difflib

# Function to read line of fasta file, returning sequence and quality scores
def seq_score(inputfile):
    fq_pair = open(inputfile, 'r')
    fq = fq_pair.read().splitlines()
    return fq[1],fq[3]
# Function to read line of template fasta
def seq(inputfile):
    with open(inputfile)as f:
        next(f)
        temp = "".join(line.strip() for line in f)
    return temp
# Function to create list of characters
def cr_item(string1):
    return [x for x in string1]

# Function to create list of index numbers for base pairs
def cr_index(r1):
    return [x for x in range(0,r1+1)]

# Function to create dictionary
def create_dict(key,val):
    return dict(zip(key,val))

# Function to translate basepairs into amino acid sequences
def translate(seq):
    table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]

    return protein
# Load file and create seuqnece score lists
os.chdir('consensus')
inputfile = ('0102.fasta')
templatefile = ('/home/longy/tam-lq.fasta')

temp_seq = seq(templatefile)
temp_protein = translate(temp_seq)
read_seq,read_sc = seq_score(inputfile)
read_protein = translate(read_seq)
basepair = cr_item(read_seq)
score = cr_item(read_sc)
num = cr_index(len(basepair))

seq_pair = create_dict(num,basepair)
score_pair = create_dict(num,score)
output_list = [li for li in difflib.ndiff(temp_protein, read_protein) if li[0] != ' ']
print(output_list)
