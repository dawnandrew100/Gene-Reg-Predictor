# Gene Regulator Predictor

Final Project for BINF 761

## Running Scripts

Python scripts can be run using `uv run <file name>.py` from the same directory
of the file.

## Pipeline Details

### Data Pre-processing

The following scripts are in the `data_processing` folder.

1. `data_filter_explore.py` opens the raw data file `mouse_Regulator_Gene.txt`
and filters the columns to only keep Regulator Symbol, Target Symbol,
and Regulator Type.
The file also filters the rows to only keep regulators that are miRNA and lncRNA.
Results are saved to `mouse_Regulator_Gene.csv`.

2. `get_sequence.py` extracts sequences from a local download of mouse genes from
NCBI `rna.fna` and a local download of mouse proteins from [STRING](https://string-db.org/cgi/organisms)
`10090.protein.sequences.v12.0.fa`. To ensure the accession number matches the
symbols in the regulator gene csv, `10090.protein.aliases.v12.0.txt` is used.
These sequences are matched with the symbols from `mouse_Regulator_Gene.csv` and
saved to `seq_info.json`. An updated CSV file containing only symbols with
accompanying sequences is created called `sequence_mouse_Regulator_Gene.csv`.

## Data Disclosure

The Gene Regulation data used in this pipeline was acquired from
[RegNetwork](http://www.zpliulab.cn/RegNetwork/download). The specific data was
the `Regulator-Gene` data for the `Complete data in mouse` dataset.

The sequence data was taken from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000001635.27/).
