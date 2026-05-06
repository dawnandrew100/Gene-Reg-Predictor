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
    # cDNA Chaos Game Representation
    radius = 10.0
    dna_points = cg.generate_square_points(4, radius=radius, rotation_deg=45)
    dna_labels = "TGCA"
    dna_labelled = cg.label_chaos_points(dna_labels, dna_points)
    dna_labelled["N"] = (0.0, 0.0)
    dna_labelled["n"] = (0.0, 0.0)
    print(dna_labelled)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # DNA sequence embeddings
    print("Opening dna json file")
    seq_file_path = "../../dna_seq_info.json"
    with open(seq_file_path) as file:
        dna_seq_dict = json.load(file)

    print("Get cDNA embeddings")
    data = [
        {"label": label, "seq": seq}
        for label, sequences in dna_seq_dict.items()
        for seq in sequences
    ]

    dna_labels = [item["label"] for item in data]
    dna_sequences = [item["seq"] for item in data]

    fcgr_list = []
    valid_dna_labels = []
    for i in trange(0, len(dna_sequences), desc="Generating DNA Tensors"):
        cgr = cg.create_chaos(dna_sequences[i], dna_labelled)
        if not cgr:
            continue
        valid_dna_labels.append(dna_labels[i])
        fcgr = cg.cgr_to_fcgr(cgr, resolution=10, radius=radius)
        fcgr_list.append(fcgr)
    full_embedding_tensor = cg.fcgr_to_tensor(fcgr_list)

    grouped_embeddings = defaultdict(list)
    for label, emb in tqdm(
        zip(valid_dna_labels, full_embedding_tensor),
        total=len(dna_labels),
        desc="Grouping DNA",
    ):
        grouped_embeddings[label].append(emb)

    # One embedding per label
    dna_final_label_vectors = {}
    for label, emb_list in tqdm(
        grouped_embeddings.items(), desc="Calculating DNA Centroids"
    ):
        # Stack the list of tensors and take the mean across the first dimension
        dna_final_label_vectors[label] = torch.stack(emb_list).mean(dim=0)

    # Save results
    torch.save(dna_final_label_vectors, "dna_embeddings.pt")

    print("Processing complete. Files saved in embeddings/")


if __name__ == "__main__":
    main()
