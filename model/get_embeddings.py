# Third party packages
import esm
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM
from tqdm import tqdm, trange

# Built ins
from collections import defaultdict
import json

"""
Reminder to try embedding with simpler custom model rather than multi-million
parameter pre-trained models
"""

def main():
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

    # Load ESM-2 model
    # 33 layer ESM-2 model with 650M params, trained on UniRef50
    model_esm, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    # alphabet.get_batch_converter returns BatchConverter class
    batch_converter = alphabet.get_batch_converter()
    model_esm.eval()  # disables dropout for deterministic results
    model_esm.to(device)

    all_prot_embeddings = []
    all_prot_labels = []
    prot_batch_size = 12
    for i in trange(0, len(data_prot), prot_batch_size, desc="ESM-2 Batches"):
        batch = data_prot[i : i + prot_batch_size]
        batch_labels, _, batch_tokens = batch_converter(batch)
        batch_tokens = batch_tokens.to(device)
        with torch.no_grad():
            results = model_esm(batch_tokens, repr_layers=[33], return_contacts=False)
            token_representations = results["representations"][33]
            batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

        # Generate per-sequence representations via averaging
        all_prot_labels.extend(batch_labels)
        for j, tokens_len in enumerate(batch_lens):
            # First and last token are beginning-of-sequence token and
            # end-of-sequence token respectively
            seq_emb = token_representations[j, 1 : tokens_len - 1].mean(dim=0).cpu()
            all_prot_embeddings.append(seq_emb)

    full_prot_tensor = torch.stack(all_prot_embeddings)
    grouped_embeddings = defaultdict(list)
    for label, emb in tqdm(
        zip(all_prot_labels, full_prot_tensor),
        total=len(all_prot_labels),
        desc="Grouping Protein Embeddings",
    ):
        grouped_embeddings[label].append(emb)

    protein_final_label_vectors = {
        label: torch.stack(emb_list).mean(dim=0)
        for label, emb_list in grouped_embeddings.items()
    }

    print(f"Embedded {len(protein_final_label_vectors)} unique protein symbols.")
    print("Protein sequence embedded")

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
    # Load GROVER model
    tokenizer = AutoTokenizer.from_pretrained("PoetschLab/GROVER")
    model_grover = AutoModelForMaskedLM.from_pretrained("PoetschLab/GROVER")
    model_grover.eval().to(device)

    all_dna_embeddings = []
    all_dna_labels = []
    dna_batch_size = 12
    for i in trange(0, len(dna_sequences), dna_batch_size, desc="GROVER Batches"):
        batch_seqs = dna_sequences[i : i + dna_batch_size]
        batch_labs = dna_labels[i : i + dna_batch_size]
        inputs = tokenizer(
            batch_seqs,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=1024,
        ).to(device)

        with torch.no_grad():
            outputs = model_grover(**inputs)
            # (Batch, Seq_Len, Hidden_Dim)
            last_hidden = outputs.last_hidden_state
            # (Batch, Seq_Len, Hidden_Dim) -> (Batch, Hidden_Dim)
            embeddings = last_hidden.mean(dim=1).cpu()

        all_dna_embeddings.append(embeddings)
        all_dna_labels.extend(batch_labs)

    full_embedding_tensor = torch.cat(all_dna_embeddings, dim=0)
    print(f"Final tensor shape: {full_embedding_tensor.shape}")

    grouped_embeddings = defaultdict(list)
    for label, emb in tqdm(
        zip(all_dna_labels, full_embedding_tensor),
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

    print(f"Embedded {len(dna_final_label_vectors)} unique DNA labels.")
    print(f"Final DNA Tensor Shape: {last_hidden.shape}")

    # Save results
    torch.save(protein_final_label_vectors, "./embeddings/protein_embeddings.pt")
    torch.save(dna_final_label_vectors, "./embeddings/dna_embeddings.pt")

    print("Processing complete. Files saved in ./embeddings/")


if __name__ == "__main__":
    main()
