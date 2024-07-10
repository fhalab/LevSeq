###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

from levseq.globals import DORADO_MODELS
import os
import glob
import subprocess


class Basecaller:

    def __init__(self, model, file_folder, save_folder, fastq = True):
        self.model = model
        self.file_folder = file_folder
        self.save_folder = save_folder
        self.fastq = fastq
        self.model_path = self.check_model()

    def check_model(self):
        """Checks if the required dorado model is downloaded."""

        model_path = os.path.join(os.path.dirname(__file__), "dorado_models")

        # Check if model already exists
        models = glob.glob(model_path + "/*")
        if DORADO_MODELS[self.model] in models:
            return os.path.join(model_path, DORADO_MODELS[self.model])
        
        else:
            self.download_model()
            model_path = os.path.join(model_path, DORADO_MODELS[self.model])

        return model_path
    
    def download_model(self):
        """Downloads the required dorado model."""
        model = DORADO_MODELS[self.model]
        
        input = f"dorado download --model {model} --directory {self.model_path}"
        subprocess.run(input, shell=True)
        return "Model downloaded."
    
    def run_dorado(self):
        """
        Runs dorado basecaller.

        Args: 
            - .pod5 files
        Returns: 
            - .bam files
        """

        model_path = self.check_model()

        if self.fastq:
            input = f"dorado basecaller {model_path} {self.file_folder} --emit-fastq > {self.save_folder}/basecalled.fastq"
        
        else:
            input = f"dorado basecaller {model_path} {self.file_folder} > {self.save_folder}/basecalled.fastq"

        subprocess.run(input, shell=True)
        
        return "Basecalling submitted"



