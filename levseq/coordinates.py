import esm
import torch
import pandas as pd
from sklearn.decomposition import PCA
import os
import argparse

def preprocess_sequence(sequence):
    """
    Preprocesses the amino acid sequence by removing everything after the first '*' (stop codon).
    """
    if '*' in sequence:
        sequence = sequence.split('*')[0]  # Take everything before the first '*'
    return sequence

def process_file(input_file, output_file=None):
    # Load the dataset
    data = pd.read_csv(input_file)

    # Remove the "Unnamed: 0" column if it exists
    if 'Unnamed: 0' in data.columns:
        data = data.drop(columns=['Unnamed: 0'])

    # Create the ID column as the combination of `Plate` and `Well`
    data['ID'] = data['Plate'] + '-' + data['Well']
    data = data[['ID'] + [col for col in data.columns if col != 'ID']]  # Reorder to make ID the first column

    # Filter valid sequences from the `aa_sequence` column
    valid_sequences = data['aa_sequence'].dropna()
    valid_sequences = valid_sequences[~valid_sequences.str.contains('#N.A.#|Deletion')]

    # Preprocess sequences to handle stop codons
    valid_sequences = valid_sequences.apply(preprocess_sequence)

    # Load the ESM-2 model
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()

    # Prepare sequences for embedding
    sequences = valid_sequences.tolist()
    sequence_names = [f"Sequence {i}" for i in range(len(sequences))]
    batch_labels, batch_strs, batch_tokens = batch_converter(list(zip(sequence_names, sequences)))

    # Extract embeddings
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33])
        embeddings = results["representations"][33]  # Use the top (last) layer representations

    # Average embeddings across residues for sequence-level representation
    sequence_embeddings = embeddings.mean(1).numpy()

    # Dimensionality Reduction using PCA
    pca = PCA(n_components=2)
    xy_coordinates = pca.fit_transform(sequence_embeddings)

    # Add x, y coordinates back to the dataframe
    xy_df = pd.DataFrame(xy_coordinates, columns=['x_coordinate', 'y_coordinate'], index=valid_sequences.index)
    data = pd.concat([data, xy_df], axis=1)

    # Determine output file location
    if output_file is None:
        input_name, input_ext = os.path.splitext(input_file)
        output_file = f"{input_name}_xy{input_ext}"

    # Save the updated dataframe to a file
    data.to_csv(output_file, index=False)
    print(f"Processed data with x, y coordinates saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate x, y coordinates for amino acid sequences")
    parser.add_argument('input_file', type=str, help="Path to the input CSV file")
    parser.add_argument('--output_file', type=str, default=None, help="Path to save the output CSV file (optional)")
    args = parser.parse_args()

    process_file(args.input_file, args.output_file)

