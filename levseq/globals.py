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

# Global parameters for minIon

# Parameter for Demultiplexing
SCORE_MATRIX = {
    ('A', 'A'): 96, ('A', 'C'): -316, ('A', 'G'): -192, ('A', 'T'): -369, ('A', 'N'): 0,
    ('C', 'A'): -316, ('C', 'C'): 100, ('C', 'G'): -352, ('C', 'T'): -295, ('C', 'N'): 0,
    ('G', 'A'): -192, ('G', 'C'): -352, ('G', 'G'): 98, ('G', 'T'): -329, ('G', 'N'): 0,
    ('T', 'A'): -369, ('T', 'C'): -295, ('T', 'G'): -329, ('T', 'T'): 100, ('T', 'N'): 0,
    ('N', 'A'): 0, ('N', 'C'): 0, ('N', 'G'): 0, ('N', 'T'): 0, ('N', 'N'): 0,
}  # Adapted from Guppy Barcoder

SW_ALIGN_PARAMS = {
    "start_gap1": 40,
    "end_gap1": 40,
    "open_gap1": 0,
    "extend_gap1": -40,
    "start_gap2": 40,
    "end_gap2": 40,
    "open_gap2": -160,
    "extend_gap2": -160,
    "min_score_barcode_front": 60.0,
    "front_window_size": 150,
    "rear_window_size": 150,
}  # Adapted from Guppy Barcoder


# Defaul target folder names from ONT
DEFAULT_TARGETS = {"Not_basecalled": ["pod5"], "Basecalled": ["fastq_pass"]}


# Codons
CODONS = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}
