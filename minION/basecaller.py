from minION.util.globals import DORADO_MODELS, MINKNOW_PATH
import os
import glob

def check_model(dorado_type):
    """Checks if the required dorado model is downloaded."""

    model_path = os.path.join(os.path.dirname(__file__), "dorado_models")

    # Check if model already exists
    models = glob.glob(model_path + "/*")
    if DORADO_MODELS[dorado_type] in models:
        return os.path.join(model_path, DORADO_MODELS[dorado_type])
   
    else:
        download_model(dorado_type, model_path)
        model_path = os.path.join(model_path, DORADO_MODELS[dorado_type])

    return model_path

def download_model(dorado_type, model_path):
    """Downloads the required dorado model."""
    model = DORADO_MODELS[dorado_type]
    
    input = f"dorado download --model {model} --directory {model_path}"
    os.system(input)
    return "Model downloaded."

def run_dorado(model, file_folder, save_folder, fastq = True):
    """Runs dorado basecaller.
    Input: .pod5 files
    Output: .bam files"""

    model_path = check_model(model)

    if fastq:
        input = f"dorado basecaller {model_path} {file_folder} --emit-fastq > {save_folder}/basecalled.fastq"
    
    else:
        input = f"dorado basecaller {model_path} {file_folder} > {save_folder}/basecalled.fastq"

    os.system(input)
    return "Basecalling submitted"




