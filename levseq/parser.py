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

import argparse
from pathlib import Path


def create_parser(): 
    parser = argparse.ArgumentParser(description='evSeq levseq pipeline. Enter the experiment name from your run')
    
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

    parser.add_argument('--output_path',
                        metavar='o',
                        default=None,
                        type=Path,
                        required=False,
                        help='Path to output folder. If not given, the output folder will be created in the current directory')

    parser.add_argument("--output_name",
                        metavar='on',
                        type=str,
                        help="Name of the output folder. If not given, the name will be the same as the experiment name")
    
    parser.add_argument("--skip_basecalling", 
                        action="store_true", 
                        help="Skip the basecalling step.")

    parser.add_argument("--skip_demultiplex", 
                        action="store_true", 
                        help="Skip the demultiplexing step.")
    
    parser.add_argument("--skip_consensus",
                        action="store_true",
                        help="Skip the consensus step.")

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