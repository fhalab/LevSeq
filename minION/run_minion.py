import os
import json
from pathlib import Path
import concurrent.futures
from minION.util.parser import create_parser, check_parser
from minION.util.IO_processor import create_folder, find_experiment_folder, find_folder, get_barcode_dict, concatenate_fastq_files, find_experiment_files
from minION.basecaller import run_dorado, check_model
from minION.demultiplexer import run_demultiplexer
from minION.consensus import process_fastq, consensus_prompt, run_medaka, medeka_stitch_prompt
from minION.analyser import get_variant_df
from minION.util.globals import BARCODES, MEDAKA_MODELS, DEFAULT_TARGETS



def process_barcode(reverse, forward, ref_seq):
    # Concat all fastq files
    fastq_file = process_fastq(forward)

    # Output directory
    output_dir = os.path.join(forward, "medaka")

    # Run Consensus
    prompt = consensus_prompt(fastq_file, output_dir, ref_seq, n_threads = 4, model = "default")
    run_medaka(prompt)

    # Run Stitch
    final_consensus = os.path.join(forward, "final_consensus.fasta")
    prompt = medeka_stitch_prompt(forward, ref_seq, final_consensus, qualities = True)

    # Align and calculate the Phred Quality Score
    run_medaka(prompt)

    
def main(args, parallel = True):
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
    
    basecall_folder = result_folder / "basecalled_filtered_un"
    
    ### ----- Basecaller ----- ###

    if not args.skip_basecalling:
        pod5_files = find_folder(experiment_folder, "pod5_pass")
        basecall_folder = os.path.join(result_folder, "basecalled_filtered")
        # Create a basecall folder if not exists
        Path(basecall_folder).mkdir(parents=True, exist_ok=True)
        run_dorado(basecall_model, pod5_files, basecall_folder, fastq = True)

    # ### ----- Demultiplex ----- ###

    if not args.skip_demultiplex:

        run_demultiplexer(result_folder, BARCODES, 60, 60, basecall_folder = basecall_folder)
        
        
    demultiplex_folder = result_folder / "demultiplex_60_un"

    ### ----- Consensus ----- ###

    barcode_dict = get_barcode_dict(demultiplex_folder)

    if not args.skip_consensus:

        if parallel:

            # Using ThreadPoolExecutor to parallelize the process
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
                futures = []
                for reverse in barcode_dict.keys():
                    for forward in barcode_dict[reverse]:
                        futures.append(executor.submit(process_barcode, reverse, forward, args.ref))

                # Collecting results (if any) and handling exceptions
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as exc:
                        print(f"Generated an exception: {exc}")

        else:
            for reverse in barcode_dict.keys():
                for forward in barcode_dict[reverse]:
                    print(f"Processing {os.path.basename(forward)}")
                    # Check if consensus file already exists
                    if os.path.exists(os.path.join(forward, "final_consensus.fasta")):
                        print(f"Consensus file in {os.path.basename(forward)} already exists")
                        continue

                    # Concat all fastq files
                    fastq_file = process_fastq(forward)

                    # Output directory
                    output_dir = os.path.join(forward, "medaka")

                    # Run Consensus
                    prompt = consensus_prompt(fastq_file, 
                                              output_dir, 
                                              args.ref, 
                                              n_threads = 2,
                                              model = "default")
                    run_medaka(prompt)

                    # Run Stitch
                    final_consensus = os.path.join(forward, "final_consensus.fasta")

                    prompt = medeka_stitch_prompt(forward, args.ref, final_consensus, qualities = True)

                    # Align and calculate the Phred Quality Score
                    run_medaka(prompt)


    # ### ----- Variant Data Frame ----- ###
    
    variant_df = get_variant_df(demultiplex_folder, args.ref, sequences=True)


    filename = f'{args.experiment_name}_{basecall_model}_variant_df_50.csv'
    filepath = os.path.join(result_folder, filename)


    variant_df.to_csv(filepath, index=False)

    return "Success!"

if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    check_parser(parser, args)
    print(main(args, parallel=False))
    