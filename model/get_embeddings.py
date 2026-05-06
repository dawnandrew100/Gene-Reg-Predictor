# Third party packages
import pandas as pd
import torch
from tqdm import tqdm, trange

# Built ins
from collections import defaultdict
import json

# Local import
import chaos_game as cg


def main():
    points = generate_square_points(4, radius=10.0, rotation_deg=45)
    dna_labels = "ACGT"
    dna_labelled = label_chaos_points(labels, v)
    print(dna_labelled)
    seq = "CGATTCGACTAGTGCACTAGTGCATGCATGCATTGCAGCGATGCACTGTGCATTCACTGATGCTAGCTGA"
    cgr = create_chaos(seq, (0, 0), labelled)
    print(cgr)
    print(cgr_to_fcgr(cgr, 10, radius=10.0))

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # Protein sequence embeddings
    print("Opening protein json file")
    seq_file_path = "../protein_seq_info.json"
    with open(seq_file_path) as file:
        protein_seq_dict = json.load(file)

    print("Convert JSON to correct input shape for ESM")
    # Convert sequences to Sequence[Tuple[str, str]]
    data_prot = [
        (symbol, seq)
        for symbol, sequences in protein_seq_dict.items()
        for seq in sequences
    ]

    # DNA sequence embeddings
    print("Opening dna json file")
    seq_file_path = "../dna_seq_info.json"
    with open(seq_file_path) as file:
        dna_seq_dict = json.load(file)

    print("Convert JSON to correct input shape for GROVER")
    data = [
        {"label": label, "seq": seq}
        for label, sequences in dna_seq_dict.items()
        for seq in sequences
    ]

    dna_labels = [item["label"] for item in data]
    dna_sequences = [item["seq"] for item in data]

    # Save results
    # torch.save(protein_final_label_vectors, "./embeddings/protein_embeddings.pt")
    # torch.save(dna_final_label_vectors, "./embeddings/dna_embeddings.pt")

    print("Processing complete. Files saved in ./embeddings/")


if __name__ == "__main__":
    main()
