import os
import json
from pathlib import Path
from minION.util.parser import create_parser, check_parser
from minION.util.IO_processor import create_folder, find_experiment_folder, find_folder, get_barcode_dict, concat_all_fastq
from minION.basecaller import run_dorado, check_model
from minION.demultiplexer import run_demultiplexer
from minION.consensus import process_fastq, consensus_prompt, run_medaka, medeka_stitch_prompt
from minION.util.globals import BARCODES, MEDAKA_MODELS

    
def main(args):

    # Arguments
    experiment_name = args.experiment_name
    output_path = args.output
    output_name = args.output_name
    ref_seq = args.ref


    result_folder = create_folder(experiment_name, output_path, output_name = output_name)

    experiment_folder = find_experiment_folder(experiment_name)

    pod5_files = find_folder(experiment_folder, "pod5")

    ### ----- Basecaller ----- ###
    basecall_folder = os.path.join(result_folder, "basecalled")

    # Create a basecall folder if not exists
    Path(basecall_folder).mkdir(parents=True, exist_ok=True)

    run_dorado("fast", pod5_files, basecall_folder, fastq = True)

    ### ----- Demultiplex ----- ###    
    run_demultiplexer(result_folder, BARCODES, 60, 50)

    demultiplex_folder = os.path.join(result_folder, "demultiplex")
    
    

    ### ----- Consensus ----- ###
    barcode_dict = get_barcode_dict(demultiplex_folder)

   
    for reverse in barcode_dict.keys():
        for forward in barcode_dict[reverse]:
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


  
    return demultiplex_folder

if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    check_parser(parser, args)
    print(main(args))
    