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
import random

import numpy as np

from levseq.variantcaller import *
from Bio import AlignIO
from sklearn.decomposition import PCA

"""
Functions for interfacing with users. 

1. Converting the final sequences to MSA
2. Validation that same sequences got the same fitness
3. Info/summaries based on what the mutations were that did good
4. Encoding etc for the sequneces
5. Generating figures.

Nonpolar (Hydrophobic)
Alanine (Ala, A) - Nonpolar, aliphatic side chain
Isoleucine (Ile, I) - Nonpolar, aliphatic side chain
Leucine (Leu, L) - Nonpolar, aliphatic side chain
Methionine (Met, M) - Nonpolar, sulfur-containing side chain
Phenylalanine (Phe, F) - Nonpolar, aromatic side chain
Proline (Pro, P) - Nonpolar, cyclic aliphatic side chain
Tryptophan (Trp, W) - Nonpolar, aromatic side chain
Valine (Val, V) - Nonpolar, aliphatic side chain
Polar, Uncharged
Asparagine (Asn, N) - Polar, amide-containing side chain
Cysteine (Cys, C) - Polar, sulfur-containing side chain
Glutamine (Gln, Q) - Polar, amide-containing side chain
Serine (Ser, S) - Polar, hydroxyl-containing side chain
Threonine (Thr, T) - Polar, hydroxyl-containing side chain
Tyrosine (Tyr, Y) - Polar, aromatic side chain with a hydroxyl group
Polar, Acidic (Negatively Charged at Physiological pH)
Aspartic Acid (Asp, D) - Acidic, carboxylate-containing side chain
Glutamic Acid (Glu, E) - Acidic, carboxylate-containing side chain
Polar, Basic (Positively Charged at Physiological pH)
Arginine (Arg, R) - Basic, contains a guanidinium group
Histidine (His, H) - Basic, contains an imidazole group
Lysine (Lys, K) - Basic, contains an amino group
Special Cases
Glycine (Gly, G) - The simplest amino acid, with a hydrogen as its side chain. It's often classified as nonpolar due to its minimal side chain, but its small size allows it to fit into both polar and nonpolar environments, making it quite versatile.
"""

# Define the standard amino acids
amino_acids = 'GAVCPLIMWFKRHSTYNQDE'
# Create a mapping from amino acids to their index
aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}
# init 0 as default value

from sciutil import SciUtil


# Just pretty printing.
u = SciUtil()


def convert_variant_df_to_seqs(variant_df, parent_seq):
    """
    Converts the variant DF to a MSA.
    """
    # Get the sequence from the DF, using the reference we'll convert this to the aa and then run a MSA on it
    seqs = [translate(parent_seq)]  # always start with the parent
    seq_ids = ['PARENT']
    for plate, well, predicted_variant in variant_df[['Plate', 'Well', 'Variant']].values:
        label = f'{plate} {well}'
        seq = list(parent_seq)
        if isinstance(predicted_variant, str):
            try:
                for mutation in predicted_variant.split('_'):
                    if 'PARENT' in mutation:
                        continue  # We skip this!
                    else:
                        # true_variant is a sequence while predicated variant is just the mutations
                        if 'DEL' not in mutation:
                            mut_pos = int(mutation[1:-1]) - 1  # A1T
                            mut = mutation[-1]
                        else:
                            mut_pos = int(mutation[1:].replace('DEL', '')) - 1
                            mut = 'DEL'
                        # This is the mutation at that point so we add it to our seq
                        if mut == 'DEL':
                            seq[mut_pos] = '-'  # We'll remove these later
                        else:
                            seq[mut_pos] = mut
            except Exception as e:
                print(e)
                print(plate, well, mut_pos, len(predicted_variant), len(parent_seq))
            seqs.append(translate(''.join(seq).replace('-', '')))  # Remove the gaps and translate the sequence
            seq_ids.append(label)
        else:
            print(label, predicted_variant)
    return seqs, seq_ids


def get_colour(aa):
    """
    Get colour for an amino acid sequence.
    """
    return 0


def make_msa(seqs, seq_ids, file_to_align='/tmp/msa.fa'):
    """
    Potentialy change this so that the MSA file has the unique time stamp so that we don't get overriding issues.
    """
    with open(file_to_align, 'w+') as fout:
        for i, seq in enumerate(seqs):
            if isinstance(seq, str) and len(seq) > 0:
                fout.write(f'>{seq_ids[i]}\n{seq}\n')

    # Now make the msa
    msa_file = f'{file_to_align.replace(".fa", "_msa.fa")}'
    # Write each one as a fasta file then run the clustal and then the tree
    os.system(f'../software/./clustal-omega-1.2.3-macosx --force -i {file_to_align} -o {msa_file} -v')
    u.dp(['Done MSA'])
    # Reading the alignment file
    alignment = AlignIO.read(msa_file, 'fasta')
    return alignment


def make_pca(encoded_sequences):
    pca = PCA(n_components=2)
    # Fit PCA on the standardized data
    return pca.fit_transform(encoded_sequences)


def one_hot_encode(sequence):
    sequence = sequence.replace('X', '-')
    sequence = sequence.replace('?', '-')

    # Initialize an array to hold the one-hot encoded sequence
    one_hot_sequence = np.zeros((len(sequence), len(amino_acids)), dtype=int)

    # Fill the array with one-hot encodings
    for i, aa in enumerate(sequence):
        if aa_to_index.get(aa) is not None:
            one_hot_sequence[i, aa_to_index.get(aa)] = 1

    flat_seq = one_hot_sequence.flatten()
    return flat_seq
