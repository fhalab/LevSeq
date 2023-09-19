# The idea is tat whenever sommeone starts a new experiment, the script will be called and search for the name of the experiment. After it finds it, it will extract
# basic information about the experiment name, the data & so on. It will then create a json file, which then can be used for rest of the scripts. The json file will
# also be useful as a log file for the experiment. The json file will be stored in the same directory as the experiment. The script will also create a directory for

import json 
import os 
from minION.util.globals import MINKNOW_PATH



def read_json(json_file):
    """Read the json file and return the dictionary
    Input: Path to json file
    Output: Dictionary"""
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data