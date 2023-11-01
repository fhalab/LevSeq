import os
import json
from pathlib import Path
import concurrent.futures
from minION.util.parser import create_parser, check_parser
from minION.util.IO_processor import create_folder, find_experiment_folder, find_folder, get_barcode_dict, concatenate_fastq_files, find_experiment_files
from minION.basecaller import run_dorado, check_model
from minION.demultiplexer import run_demultiplexer
from minION.consensus import process_fastq,get_consensus
from minION.analyser import get_variant_df_nn
from minION.util.globals import BARCODES, MEDAKA_MODELS, DEFAULT_TARGETS


  
def main(args, parallel = False):
    """Main Function to run evSeq minION. In order to run evSeq-minION, make sure that you have run the sequencing successfully. Ideally, also check the 
    quality report from Nanopore to make sure that the sequencing was successful and the quality of the reads are good. However, this script will also provide a quality report at the end of the run. \n
    Please note that this script is designed to run on a Linux machine. 
    
    Args:
        - args (argparse.Namespace): Arguments from the command line
            - experiment_name (str): Name of the experiment. This will be used to create a folder to store the results.
            - output_path (str): Path to the output folder. Defaults to the current working directory.
            - TOML config file (str): Path to the TOML config file. Defaults to the config.toml file in the current working directory. #TODO: Add a default config file
            - output_name (str): Name of the output folder. Defaults to the experiment_name.
            - ref (str): Path to the reference sequence fasta file.
            - skip_basecalling (bool): If True, the script will skip the basecalling step. Defaults to False.
            - skip_demultiplex (bool): If True, the script will skip the demultiplexing step. Defaults to False.
            - skip_consensus (bool): If True, the script will skip the consensus step. Defaults to False.
        
        - parallel (bool, optional): If True, the script will run in parallel. Defaults to True.
        
    """

    # Arguments
    
    basecall_model = "sup"

    result_folder = create_folder(args.experiment_name, basecall_model, target_path=args.output_path, output_name=args.output_name)

    experiment_folder = find_experiment_folder(args.experiment_name)

    #TODO: Find Files based if it was basecalled or not
    for key, value in DEFAULT_TARGETS.items():
        file_path = find_experiment_files(experiment_folder, value)
    
    basecall_folder = result_folder / "basecalled_filtered"
    
    ### ----- Basecaller ----- ###

    if not args.skip_basecalling:
        pod5_files = find_folder(experiment_folder, "pod5_pass")
        basecall_folder = os.path.join(result_folder, "basecalled_filtered")
        # Create a basecall folder if not exists
        Path(basecall_folder).mkdir(parents=True, exist_ok=True)
        run_dorado(basecall_model, pod5_files, basecall_folder, fastq = True)

    # ### ----- Demultiplex ----- ###

    if not args.skip_demultiplex:

        run_demultiplexer(result_folder, BARCODES, 45, 45, basecall_folder = basecall_folder)
        
        
    demultiplex_folder = result_folder / "demultiplex_45"

    consensus_folder_name = "consensus"

    ### ----- Consensus ----- ###

    barcode_dict = get_barcode_dict(demultiplex_folder)

    if not args.skip_consensus:

        if parallel:

            #TODO : Add parallel processing
            pass

        else:
            for reverse in barcode_dict.keys():
                for forward in barcode_dict[reverse]:

                    print(f"Processing {os.path.basename(forward)}")

                    # Check if consensus file already exists
                    if os.path.exists(os.path.join(forward, consensus_folder_name , "consensus.fastq")):
                        print(f"Consensus file in {os.path.basename(forward)} already exists")
                        continue
                    
                    print("Processing fastq files")
                    concatenate_fastq_files(Path(forward), filename = "concated", prefix = "demultiplexed", delete = True)
                    print("Processing fastq files done")
    
                    # Run Consensus
                    get_consensus(Path(forward), args.ref, "consensus.fastq", qualities=True, consensus_folder=consensus_folder_name)              



    # ### ----- Variant Data Frame ----- ###
    
    variant_df = get_variant_df_nn(demultiplex_folder, args.ref, consensus_folder_name="consensus" ,sequences=True)


    filename = f'{args.experiment_name}_{basecall_model}_guppy_45.csv'
    filepath = os.path.join(result_folder, filename)


    variant_df.to_csv(filepath, index=False)

    return "Success!"

if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    check_parser(parser, args)
    print(main(args, parallel=False))
    