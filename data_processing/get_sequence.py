# Third party packages
from Bio import Entrez
import pandas as pd
from tqdm import tqdm

# Packages made by me
from biobase.parser import FastaFileParser, fasta_parser, FastaRecord

# Built ins
from collections import defaultdict
import getpass
import json
import re


def main():
    # Load all mouse data
    mouse_file_name = "mouse_Regulator_Gene"
    mouse_file_path = f"../{mouse_file_name}.csv"
    mouse_df = pd.read_csv(mouse_file_path)

    mouse_df["Regulator Symbol"] = mouse_df["Regulator Symbol"].str.replace(
        r"mmu-mi[rR]-([^-]+).*", r"Mir\1", regex=True
    )
    mouse_df["Regulator Symbol"] = mouse_df["Regulator Symbol"].str.replace(
        r"mmu-let-([^-]+).*", r"Mirlet\1", regex=True
    )
    # Initial row count: 1,881,858
    print(mouse_df)

    regulators = mouse_df["Regulator Symbol"].unique().tolist()
    targets = mouse_df["Target Symbol"].unique().tolist()

    symbols = set(regulators + targets)
    print(f"There are {len(symbols)} symbols.")
    mouse_seq_dict: dict[str, list[str]] = defaultdict(list)

    # Fetch and store records from local ncbi
    rna_file_path = "../rna.fna"
    mouse_fasta_records = iter(FastaFileParser(rna_file_path))
    mouse_seq_dict = fetch_sequences_from_local_ncbi(
        mouse_fasta_records, mouse_seq_dict, symbols
    )

    # Fetch and store remaining records from local copy of strings library
    strings_file_path = "../10090.protein.aliases.v12.0.txt"
    strings_df = pd.read_csv(strings_file_path, sep="\t")
    strings_seq_file_path = "../10090.protein.sequences.v12.0.fa"
    strings_fasta_records = iter(FastaFileParser(strings_seq_file_path))
    mouse_seq_dict = fetch_sequences_from_local_strings(
        strings_df, strings_fasta_records, mouse_seq_dict, symbols
    )

    ## I tried querying ncbi for the remaining sequences but ncbi took an incredibly long time to fetch sequences
    # remaining_recs = symbols - set(mouse_seq_dict.keys())
    # mouse_seq_dict = fetch_sequences_from_ncbi(list(remaining_recs), mouse_seq_dict)

    recs_with_seqs = set(mouse_seq_dict.keys())
    filtered_df = mouse_df[
        mouse_df["Regulator Symbol"].isin(recs_with_seqs)
        & mouse_df["Target Symbol"].isin(recs_with_seqs)
    ]
    # Final row count: 363,671
    print(filtered_df)
    print(f"Finished with {len(recs_with_seqs):,} symbols")
    # Started with 29,745 unique symbols
    # 22,138 symbols after local ncbi sequences
    # 25,869 symbols after local strings sequences
    # Over a million rows were removed from the df 1,881,859 -> 1,368,611
    print(filtered_df.describe())
    print(filtered_df["Regulator Type"].value_counts())

    # Write filtered data to new CSV file
    filtered_filepath = f"../sequence_{mouse_file_name}.csv"
    filtered_df.to_csv(filtered_filepath, index=False)

    sequence_info_file = "../seq_info.json"
    with open(sequence_info_file, "w") as file:
        json.dump(mouse_seq_dict, file, indent=4)


def extract_symbol(record: FastaRecord) -> str | None:
    if match := re.search(r"(?<=\()([^-\)]+)", record.name):
        return match.group()
    return None


def fetch_sequences_from_local_ncbi(
    fasta: Iterator[FastaRecord], record_dict: dict[str, list[str]], symbols: set[str]
) -> dict[str, list[str]]:
    print("Extracting sequences")
    # Takes about 1-2 seconds to run without translation
    # Processes a total of 136,225 records at about 70k per second
    for record in tqdm(fasta, desc="Processing FASTA", unit=" records"):
        record.name = extract_symbol(record)
        if record.name in symbols:
            record_dict[record.name].append(record.seq)

    # 100,728 sequences for 22,614 records
    print(f"Finished. Found sequences for {len(record_dict):,} records.")

    return record_dict


def fetch_sequences_from_local_strings(
    df: pd.DataFrame,
    fasta: Iterator[FastaRecord],
    record_dict: dict[str, list[str]],
    symbols: set[str],
) -> dict[str, list[str]]:

    filtered_df = df[
        df["alias"].isin(symbols) & ~df["alias"].isin(set(record_dict.keys()))
    ]
    unique_df = filtered_df.drop_duplicates(subset="#string_protein_id", keep="first")
    print(df)
    print(filtered_df)
    print(unique_df)
    id_to_symbol_dict = dict(zip(unique_df["#string_protein_id"], unique_df["alias"]))
    # This happens in less than a second thanks to dictionary map
    for record in tqdm(fasta, desc="Processing FASTA", unit=" records"):
        symbol = id_to_symbol_dict.get(record.id)
        if symbol:
            record_dict[symbol].append(record.seq)
    return record_dict


def fetch_sequences_from_ncbi(
    symbols: list[str], records: dict[str, list[str]], batch_size: int = 50
):
    Entrez.email = input("Enter your email: ").strip()
    Entrez.api_key = getpass.getpass(
        "Enter your NCBI API Key (press Enter to skip): "
    ).strip()

    for i in tqdm(range(0, len(symbols), batch_size), desc="Downloading from NCBI"):
        batch = symbols[i : i + batch_size]
        query = " OR ".join([f"{sym}[Gene Name] AND mouse[Organism]" for sym in batch])

        try:
            # Search for the IDs
            search_handle = Entrez.esearch(
                db="nucleotide", term=query, retmax=batch_size
            )
            search_results = Entrez.read(search_handle)
            id_list = search_results.get("IdList", [])
            search_handle.close()

            if not id_list:
                continue

            # Fetch FASTA sequences
            fetch_handle = Entrez.efetch(
                db="nucleotide", id=id_list, rettype="fasta", retmode="text"
            )
            data = fetch_handle.read()
            fetch_handle.close()

            # Parse text data into records
            recs = fasta_parser(data)
            for rec in recs:
                records[rec.name].append(rec.seq)

        except Exception as e:
            print(f"Error fetching batch starting at {i}: {e}")
    return records


if __name__ == "__main__":
    main()
