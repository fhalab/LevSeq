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
"""
Contain argument parsers used for command line interface and web interface
"""
# Import packages
import os
import tqdm
import argparse
import pandas as pd

# Import local packages
from levseq.run_levseq import run_LevSeq

# Get the working directory
CWD = os.getcwd()

# Set default arguments 
padding_start = 0
padding_end = 0
min_depth = 5
threshold = 0.2
basecall_model = 'sup'


# Build the CLI argparser
def build_cli_parser():
    # Initialize
    parser = argparse.ArgumentParser()

    # Add required arguments
    required_args_group = parser.add_argument_group("Required Arguments", "Arguments required for each run")
    required_args_group.add_argument('name',
            help = 'User defined name for the output folder')
    required_args_group.add_argument("path",
            help="Path to folder containing fastq.pass or pod5_pass files.")
    required_args_group.add_argument("summary",
            help="CSV file containig barcodes used, name of each plate and reference sequence in string")
    # Add optional arguments
    optional_args_group = parser.add_argument_group("Optional Arguments", "Aditional arguments")
    optional_args_group.add_argument("--output",
                                     help="Save location for run. Defaults to current working directory.",
                                     required=False,
                                     default=CWD)
    optional_args_group.add_argument("--perform_basecalling",
                                     action="store_true",
                                     help="Skip the basecalling step, default is false")
    optional_args_group.add_argument("--skip_demultiplexing",
                                     action="store_true",
                                     help="Skip the demultiplexing step, default is false")
    optional_args_group.add_argument("--skip_variantcalling",
                                     action="store_true",
                                     help="Skip the variant calling step, default is false")
    optional_args_group.add_argument("--oligopool",
                                     action="store_true",
                                     help="Whether this experiment came from an oligopool, default is false.")
    optional_args_group.add_argument("--show_msa",
                                     default=False,
                                     help="Skip showing msa")  
    # if cl_args.get('fitness_files') and cl_args.get('smiles'):
    optional_args_group.add_argument("--fitness_files",
                                    default=None,
                                    help="A comma separated list of fitness files (full path) with string quotation marks around them.")
    optional_args_group.add_argument("--smiles",
                                default=None,
                                help="A smiles string of the reaction with quotation marks around.")
    optional_args_group.add_argument("--compound",
                            default=None,
                            help="The compound in the fitness files (e.g. pDT or pdt - case sensitive).")     
    optional_args_group.add_argument("--variant_df",
                        default=None,
                        help="The variant dataframe to combine with fitness data.")                                   
    return parser


def combine_seq_func_data(cl_args):
    # Also check if we have any fitness data 
    if cl_args.get('fitness_files') and cl_args.get('smiles') and cl_args.get('variant_df'):
        variant_filename = cl_args.get('variant_df')
        variant_df = pd.read_csv(variant_filename)
        # Combine the fitness data with the plate data (note the barcode has to be the last _[barcode])
        # The smiles has to be the reaction smiles
        function_files = cl_args.get('fitness_files')
        compound_name = cl_args.get('compound') if cl_args.get('compound') else 'pdt'
        print(function_files, compound_name)
        all_function_df = pd.DataFrame()
        for function_file in function_files.split(','):
            barcode = function_file.split('.csv')[0].split('_')[-1]
            function_df = pd.read_csv(f'{function_file}')
            function_df.columns = [c.replace('\n', ' ') for c in function_df.columns]
            function_df['function_well'] = [x.split('-')[-1] if isinstance(x, str) else None for x in function_df['Sample Vial Number'].values]
            function_df['function_barcode_plate'] = barcode
            function_df = function_df[function_df['Compound Name'] == compound_name] # We only use pdt or Pdt
            # Convert it to numeric 
            function_df['Area'] = pd.to_numeric(function_df['Area'], errors='coerce')

            function_df['barcode_well'] = [f'{p}_{w}' for w, p in function_df[['function_well', 'function_barcode_plate']].values]
            function_df['filename'] = function_file
            print(function_df.head())
            all_function_df = pd.concat([all_function_df, function_df])
        # Join this with the variant_df barcode plate
        variant_df['barcode_well'] = [f'{p}_{w}' for w, p in variant_df[['Well', 'barcode_plate']].values]
        # Join the two
        variant_df.set_index('barcode_well', inplace=True)
        all_function_df.set_index('barcode_well', inplace=True)
        variant_df = variant_df.join(all_function_df, how='left')
        reaction_smiles = cl_args.get('smiles')
        variant_df['smiles_string'] = reaction_smiles.split('>>')[-1]
        variant_df['reaction_smiles'] = reaction_smiles
        variant_df.columns = [c.lower().replace(' ', '_') for c in variant_df.columns]
        variant_df.rename(columns={'area': 'fitness_value'}, inplace=True)
        variant_df.to_csv(f'{variant_filename.replace(".csv", "_seqfunc.csv")}')
        
        # levseq levseq_4.1 ref.csv fitness --fitness_files "20250712_epPCR_Q06714_37.csv,20250712_epPCR_Q06714_38.csv,20250712_epPCR_Q06714_39.csv,20250712_epPCR_Q06714_40.csv" --smiles 'O=P(OC1=CC=CC=C1)(OC2=CC=CC=C2)OC3=CC=CC=C3>>O=P(O)(OC4=CC=CC=C4)OC5=CC=CC=C5'  --compound dPPi --variant_df visualization_partial.csv
        return variant_df
        
# Execute LevSeq
def execute_LevSeq():
    # Build parser
    parser = build_cli_parser()
    # Parse the arguments
    CL_ARGS = vars(parser.parse_args())
    if CL_ARGS.get('fitness_files') and CL_ARGS.get('smiles') and CL_ARGS.get('variant_df'):
        print('Combining LevSeq')
        return combine_seq_func_data(CL_ARGS)
    # Set up progres bar
    tqdm_fn = tqdm.tqdm
    # Run LevSeq
    try:
        from levseq import __version__
        print(f"Starting LevSeq v{__version__}...")
        run_LevSeq(CL_ARGS, tqdm_fn)
        print(f"Run completed successfully. Results and logs stored in {os.path.join(CL_ARGS.get('output', CWD), CL_ARGS.get('name', ''))}")
    except Exception as e:
        print(f"Error: {e}")
        print(f"Check error logs for details in {os.path.join(CL_ARGS.get('output', CWD), CL_ARGS.get('name', ''))}")
