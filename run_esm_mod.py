import numpy as np
import pandas as pd
import torch
import esm
from Bio import SeqIO
import csv

input_ls = ['input_duf3494a.fasta'] # input filename(s)

output_ls = ['output_name.csv']

# Loop over the input/output files
for i in range(len(input_ls)):
    input_name = input_ls[i]
    output_name = output_ls[i]
    records = list(SeqIO.parse(input_name, "fasta"))

    # Load the ESM model and alphabet
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # Disables dropout for deterministic results

    # Process data in batches
    batch_size = 50  # Adjust this value based on your available memory
    num_batches = (len(records) + batch_size - 1) // batch_size

    sequence_representations = []  # To store representations for each sequence
    all_batch_data = []  # To store data for all batches

    for batch_num in range(num_batches):
        start_idx = batch_num * batch_size
        end_idx = min((batch_num + 1) * batch_size, len(records))

        # Prepare data for the current batch
        batch_data = [(record.id, str(record.seq)) for record in records[start_idx:end_idx]]
        all_batch_data.extend(batch_data)  # Collect data for all batches
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        # Extract per-residue representations (on CPU)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33]

        # Generate per-sequence representations via averaging
        # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
        for i, tokens_len in enumerate(batch_lens):
            sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))

        print(f"Processed batch {batch_num + 1}/{num_batches}")

    # Write data to CSV
    output_file = output_name

    with open(output_file, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # Write header
        header = [f"feature_{i}" for i in range(len(sequence_representations[0]))]
        writer.writerow(["sequence_id"] + header)

        # Write data for all batches
        for i, seq_repr in enumerate(sequence_representations):
            sequence_id = all_batch_data[i][0]
            writer.writerow([sequence_id] + seq_repr.tolist())

    print(f"Sequence representations written to {output_file}")
