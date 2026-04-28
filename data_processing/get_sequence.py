# Third party packages
from Bio import Entrez
import pandas as pd
from tqdm import tqdm

# Packages made by me
from biobase.parser import FastaFileParser, fasta_parser

# Built ins
from collections import defaultdict
import getpass
import io
import json
import re


def main():
    file_name = "filtered_mouse_Regulator_Gene"
    file_path = f"../{file_name}.csv"
    df = pd.read_csv(file_path)
    # Initial row count: 1,881,858
    print(df)

    regulators = df["Regulator Symbol"].unique().tolist()
    targets = df["Target Symbol"].unique().tolist()

    # Set of all of the regulator and target symbols
    symbols = set(regulators + targets)
    # Friendly reminder that processing increases from ~6k/sec to ~73k per sec
    # just by using a set instead of a list; ie from ~23 seconds to <2 seconds
    print(f"There are {len(symbols)} symbols.")

    mouse_records: dict[str, list[str]] = defaultdict(list)

    rna_file_path = "../rna.fna"
    parsed = iter(FastaFileParser(rna_file_path))
    print("Extracting sequences")

    # Takes about 1-2 seconds to run
    # Processes a total of 136,225 records at about 70k per second
    for record in tqdm(parsed, desc="Processing FASTA", unit=" records"):
        record.name = extract_symbol(record)
        if record.name in symbols:
            mouse_records[record.name].append(record.seq)

    # 100,728 sequences for 22,614 records
    print(f"Finished. Found sequences for {len(mouse_records):,} records.")

    # I tried querying ncbi for the remaining sequences but there were so many
    # that I got rate limited and my script froze

    # remaining_recs = symbols - set(mouse_records.keys())
    # mouse_records = fetch_sequences_from_ncbi(list(remaining_recs), mouse_records)

    recs_with_seqs = set(mouse_records.keys())
    filtered_df = df[
        df["Regulator Symbol"].isin(recs_with_seqs)
        & df["Target Symbol"].isin(recs_with_seqs)
    ]
    # Final row count: 324,566
    print(filtered_df)
    # Started with 31,703 unique symbols and ended with 22,614
    # Over a million rows were removed from the df 1,881,859 -> 324,566

    # Write filtered data to new CSV file
    filtered_filepath = f"../sequence_{file_name}.csv"
    filtered_df.to_csv(filtered_filepath, index=False)

    sequence_info_file = "../seq_info.json"
    with open(sequence_info_file, "w") as file:
        json.dump(mouse_records, file, indent=4)


def extract_symbol(record: FastaRecord) -> str | None:
    if match := re.search(r"(?<=\()\w+(?=\))", record.name):
        return match.group()
    return None


def fetch_sequences_from_ncbi(
    symbols: list[str], records: dict[str, list[str]], batch_size: int = 300
):
    Entrez.email = input("Enter your email: ").strip()
    Entrez.api_key = getpass.getpass("Enter your NCBI API Key (press Enter to skip): ").strip()

    for i in tqdm(range(0, len(symbols), batch_size), desc="Downloading from NCBI"):
        batch = symbols[i : i + batch_size]
        query = " OR ".join([f"{sym}[Gene Name] AND mouse[Organism]" for sym in batch])

        try:
            # Search for the IDs
            search_handle = Entrez.esearch(
                db="nucleotide", term=query, retmax=batch_size
            )
            search_results = Entrez.read(search_handle)
            id_list = search_results["IdList"]
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
