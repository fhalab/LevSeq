# Import

from minION.util import IO_processor
from minION import analyser
from minION import consensus

import importlib

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import gzip


# Functions 

# Pandas releated

def rewrite_Plate_name(entry):
    entry = entry[-1]

    return int(entry)

def rewrite_Well_name(entry):
    num = int(entry[1:])
    new_entry = entry[0] + str(num)

    return new_entry

def single_plate_annotation(entry):
    row = ["A", "B", "C", "D", "E", "F", "G", "H"]
    new_well_name = row[int(entry["Plate"]) - 1] + entry["Well"][1:]
    entry["Well"] = new_well_name
    return entry

def variant_combo(entry):
    if entry in ["#PARENT#", "NC"]:
        return entry
    elif entry == float("nan"):
        return entry
    else:
        new_entry = "_".join(entry)       
        return new_entry

def transform(entry):
    # If the entry is either DEAD or PARENT, return empty lists
    if entry in ["#DEAD#", "#PARENT#"]:
        return [], [], []
    else:
        # Split the entry based on the underscore
        mutations = entry.split("_")
        parent_combo = []
        positions = []
        new_aa = []
        for mutation in mutations:
            # Append the first character to parent_combo, the middle to positions, and the last to new_aa

            #Skip if parent is equal to new_aa
            if mutation[0] == mutation[-1]:
                continue

            elif mutation.find("DEL") != -1:
                parent_combo.append(mutation[0])
                positions.append(mutation[1:-3])
                new_aa.append("DEL")
            
            else:
                parent_combo.append(mutation[0])
                positions.append(mutation[1:-1])
                new_aa.append(mutation[-1])

        return parent_combo, positions, new_aa

def transform_ref(entry):

    parent_combo = []
    positions = []
    new_aa = []

    try:

        for variant in entry:

            if variant in ["#DEAD#", "#PARENT#"]:
                return ["-"], ["-"], ["-"]

            elif variant.find("DEL") != -1:
                parent_combo.append(variant[0])
                positions.append(variant[1:-3])
                new_aa.append("DEL")
            
            else:
                parent_combo.append(variant[0])
                positions.append(variant[1:-1])
                new_aa.append(variant[-1])
        return parent_combo, positions, new_aa
    
    except:
        return ["NA"], ["NA"], ["NA"]
        
def substract_by_index(positions, index):
    new_positions = []
    for pos in positions:
        new_positions.append(str(int(pos) - index))
    
    return new_positions

# Split Variant Combo
def split_variant_combo(entry):
    if entry in ["#DEAD#", "#PARENT#"]:
        return [entry]
    else:
        return entry.split("_")
    
def combine_variant_combo(entry):
    if entry in ["#DEAD#", "#PARENT#"]:
        return entry
    else:
        return "_".join(entry)

def compare_mutations(row, mutation_only = True):

    if isinstance(row['Variant_x'], float):
        return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

    if mutation_only and row['Variant_Combo_Ref'] in [["#DEAD#"], ["#PARENT#"]]:
        return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

    
    mutations1 = set(row['Variant_Combo_Ref'])
    mutations2 = set(row['Variant_x'])
    
    correct_mutations = mutations1.intersection(mutations2)
    missed_mutations = mutations1.difference(mutations2)
    
    return pd.Series([len(correct_mutations), len(missed_mutations), list(correct_mutations), list(missed_mutations)], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

def compare_only_Parent(row):

    if isinstance(row['Variant_x'], float):
        return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])
    
    if row['Variant_Combo_Ref'] in [["#PARENT#"]]:
        mutations1 = set(row['Variant_Combo_Ref'])
        mutations2 = set(row['Variant_x'])
        
        correct_mutations = mutations1.intersection(mutations2)
        missed_mutations = mutations1.difference(mutations2)
        return pd.Series([len(correct_mutations), len(missed_mutations), list(correct_mutations), list(missed_mutations)], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

    else:
         return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

def num_mutations(row):
    if isinstance(row['Variant_x'], float):
        return 0
    
    if row['Variant_Combo_Ref'] in [["#DEAD#"], ["#PARENT#"]]:
        return 0
    
    else:
        return len(row['Variant_Combo_Ref'])
    
def count_parent(row):
    if isinstance(row['Variant_x'], float):
        return 0
    
    if row['Variant_Combo_Ref'] in [["#PARENT#"]]:
        return 1
    
    else:
        return 0

def count_all(row):

    if isinstance(row['Variant_Combo_Ref'], float):
        return 0
    
    elif row['Variant_Combo_ref'] in [["#PARENT#"]]:
        return 1
    
    else:
        return len(row['Variant_Combo_Ref'])


# Function to create VariantCombo
def create_variant_combo(row):
    variant_combos = []
    for parent, pos, aa in zip(row['ParentCombo'], row['Positions'], row['NewAA']):
        if parent and pos and aa:
            variant_combos.append(f"{parent}{pos}{aa}")
    return variant_combos if variant_combos else ["#PARENT#"] if not row['ParentCombo'] else ["#DEAD#"]
    




