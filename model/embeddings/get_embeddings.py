# Third party packages
import pandas as pd
import torch
from tqdm import tqdm, trange

# Built ins
from collections import defaultdict
import json

# Local import
import chaos_game as cg
from embedding_nn import NeuralBedsCNN


def main():
    # cDNA Chaos Game Representation
    radius = 10.0
    res = 10
    dna_points = cg.generate_square_points(4, radius=radius, rotation_deg=45)
    dna_labels = "TGCA"
    dna_labelled = cg.label_chaos_points(dna_labels, dna_points)
    dna_labelled.update({"N": (0.0, 0.0), "n": (0.0, 0.0)})
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
    for label, sequences in tqdm(
        dna_seq_dict.items(), total=len(dna_seq_dict), desc="Generating DNA Tensors"
    ):
        for seq in sequences:
            cgr = cg.create_chaos(seq, dna_labelled)
            if cgr:
                valid_dna_labels.append(label)
                fcgr_list.append(cg.cgr_to_fcgr(cgr, resolution=res, radius=radius))

    # Convert to single [N, 1, 10, 10] tensor
    fcgr_tensor = cg.fcgr_to_tensor(fcgr_list).to(device).float()

    model = NeuralBedsCNN(embedding_dim=512).to(device)
    model.eval()
    all_embeddings = []
    batch_size = 64

    with torch.no_grad():
        for i in trange(0, len(fcgr_tensor), batch_size, desc="Inferring Embeddings"):
            batch = fcgr_tensor[i : i + batch_size].float().to(device)
            emb = model(batch, return_embedding=True)
            all_embeddings.append(emb.cpu())

        all_embeddings = torch.cat(all_embeddings, dim=0)

    grouped_embeddings = defaultdict(list)
    for label, emb in tqdm(
        zip(valid_dna_labels, all_embeddings),
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
