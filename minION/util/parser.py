import argparse
from pathlib import Path


def create_parser(): 
    parser = argparse.ArgumentParser(description='evSeq minION pipeline. Enter the experiment name from your run')
    
    # Arguments
    parser.add_argument('--experiment_name',
                        metavar='n', 
                        type=str,
                        required=True,
                        help='Name of experiment. The name must overlap with the name given for Sequencing')
    
    parser.add_argument('--ref',
                        metavar='r',
                        type=Path,
                        required=True,
                        help='Path to reference sequence.')

    parser.add_argument('--output',
                        metavar='o',
                        type=Path,
                        required=False,
                        help='Path to output folder. If not given, the output folder will be created in the current directory')

    parser.add_argument("--output_name",
                        metavar='on',
                        type=str,
                        help="Name of the output folder. If not given, the name will be the same as the experiment name")

    parser.add_argument('--json_file',
                        metavar='j',
                        type=Path,
                        required=False,
                        help='Path to json file. If not given, default json file will be used')

    
    return parser


def check_parser(parser, args):
    """Check if the parser argumemtns are valid
    Input: parser object, parser arguments
    Output: True if valid, else exit"""

    if args.experiment_name is None:
        parser.print_help()
        assert "Please enter the experiment name"
        exit(1)
    else:
        return True