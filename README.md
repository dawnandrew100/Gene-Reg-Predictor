# Gene Regulator Predictor

Final Project for BINF 761

## Running Scripts

Python scripts can be run using `uv run <file name>.py` from the same directory
of the file.

## Pipeline Details

### Data Pre-processing

The following scripts are in the `data_pre-processing/` folder.

1. `data_filter_explore.py` opens the raw data file `mouse_Regulator_Gene.txt`
and filters the columns to only keep Regulator Symbol, Target Symbol,
and Regulator Type.
The file also filters the rows to only keep regulators that are miRNA and lncRNA.

- Input files
  - `mouse_Regulator_Gene.txt`
- Ouput files
  - `mouse_Regulator_Gene.csv`

2. `get_sequence.py` extracts sequences from a local download of mouse genes from
NCBI `rna.fna` and a local download of mouse proteins from [STRING](https://string-db.org/cgi/download?species=Mus+musculus)
`10090.protein.sequences.v12.0.fa`. To ensure the accession number matches the
symbols in the regulator gene csv, `10090.protein.aliases.v12.0.txt` is used.
These sequences are matched with the symbols from `mouse_Regulator_Gene.csv` and
saved to `dna_seq_info.json` and `protein_seq_info.json`. An updated CSV file
containing only symbols with accompanying sequences is created called
`sequence_mouse_Regulator_Gene.csv`.

- Input files
  - `mouse_Regulator_Gene.csv`
  - `rna.fna`
  - `10090.protein.sequences.v12.0.fa`
  - `10090.protein.aliases.v12.0.txt`
- Ouput files
  - `dna_seq_info.json`
  - `protein_seq_info.json`
  - `sequence_mouse_Regulator_Gene.csv`
- Model
  - [Convolutional NN embedding model](https://github.com/dawnandrew100/Gene-Reg-Predictor/blob/main/model/embeddings/embedding_nn.py)
    - Generates higher dimensional embeddings of input sequences based on frequency chaos game representation (FCGR)
    - ```py
      class NeuralBedsCNN(nn.Module)
      ```

The following script is in the `model/embeddings/` folder

3. `get_embeddings.py` opens the previously generated JSON files
and converts the sequences into their respective embeddings.
The sequences are converted to a [Chaos Game](https://en.wikipedia.org/wiki/Chaos_game)
Representation (CGR) and this CGR is converted to a Frequency CGR (FCGR) by
counting the number of dots in each grid section. This FCGR represents a k-mer
size of 3.32 k-mers given by the equation $k = \frac{\log_{2} (q)}{2}$ where
k is the k-mer size and q is the number of quadrants.
The benefit of the FGCR is that sequences of all sizes are represented by a fixed
size numpy 2D array. This eliminated the need for padding. The Chaos Game functions are
located in `chaos_game.py`. The model used to produce the embeddings is is `embeddings_nn.py`.

- Input files
  - `dna_seq_info.json`
- Ouput files
  - `dna_embeddings.pt`

> [!NOTE]
> `get_embeddings_pretrained.py` is an alternative to `get_embeddings.py` but due
> to the size of the models, each model would take 4-5 days to run. The FCGR alternative
> takes about 8 minutes to run.

The following script is in the `model/` folder

4. `nn.py` opens sequence information from `sequence_mouse_Regulator_Gene.csv` and loads embeddings from
`dna_embeddings.pt` then converts embeddings into an appropriately shaped dataset using the `get_dataset` function.
This function combines regulators and targets based on the `.csv` to create positive samples then generates negative
samples on the fly by combing targets with regulators not explicitly paired in the `.csv`. The data is then split into
training and testing data and the relevant model is trained to 100 epochs.
- Model
  - [Recurrent NN predictor model](https://github.com/dawnandrew100/Gene-Reg-Predictor/blob/main/model/nn_model.py)
    - Predicts whether RNA `x` regulates gene `Y` and generates a probability of this regulatory pair
    - ```py
      class GeneRegRNNModel(nn.Module)
      ``` 
  - [Transformer predictor model](https://github.com/dawnandrew100/Gene-Reg-Predictor/blob/main/model/nn_model.py)
    - Predicts whether RNA `x` regulates gene `Y` and generates a probability of this regulatory pair
    - ```py
      class GeneRegTransformerModel(nn.Module)
      ```

## Results

The performance of the Gated Recurrent Unit (GRU) RNN and Transformer models was evaluated on gene regulation
prediction capability across four metrics; accuracy, precision, recall, and AUC-ROC.
The Transformer architecture consistently outperformed the GRU RNN across all metrics.

![RNN Loss Curve](https://github.com/dawnandrew100/Gene-Reg-Predictor/blob/main/figures/RNN_Loss_Curve.png)

![Transformer Loss Curve](https://github.com/dawnandrew100/Gene-Reg-Predictor/blob/main/figures/Transformer_Loss_Curve.png)

| Model            | Dataset  | Accuracy | Precision | Recall | AUC-ROC |
| ---------------- | -------- | -------- | --------- | ------ | ------- |
| **GRU RNN**      | Training | 0.7007   | 0.7212    | 0.7361 | 0.7658  |
|                  | Testing  | 0.6993   | 0.7198    | 0.7351 | 0.7645  |
| **Transformers** | Training | 0.7236   | 0.7357    | 0.7705 | 0.7990  |
|                  | Testing  | 0.7242   | 0.7361    | 0.7713 | 0.7995  |

* **GRU RNN Model Analysis:** The loss curve leveled off early (starting at epoch 40 and plateauing completely by epoch 60),
indicating early convergence but ultimately somewhat underfitting the data as indicated by the better than random yet still
low scores of around 70%. Even though this was the case, the model still generalised well as there were no gaps between
training and testing metrics. This suggests that while the GRU captured broad regulator-target interactions, the recurrent
architectures may have struggled with the long-distance sequence matches necessary to capture relevant motifs.
  
* **Transformer Model Analysis:** The Transformer's loss curve did not plateau within 100 epochs, indicating further room
for optimisation (e.g., via a higher learning rate). Although the metrics were overall better with the transformer model,
they could still be improved to all be >85% through hyperparameter tuning. Despite this, it achieved a substantial jump
in recall and AUC-ROC.

## Data Disclosure

The Gene Regulation data used in this pipeline was acquired from
[RegNetwork](http://www.zpliulab.cn/RegNetwork/download). The specific data was
the `Regulator-Gene` data for the `Complete data in mouse` dataset.

The sequence data was taken from
[NCBI](https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000001635.27/) for RNA/DNA sequences and
[STRING](https://string-db.org/cgi/download?species=Mus+musculus) for protein sequences.

## Sequence Embeddings

- Pre-trained models
  - Protein sequence embedding
    - [Meta's Evolutionary Scale Modeling](https://github.com/facebookresearch/esm)
    (esm) model
    - [Model article](https://www.biorxiv.org/content/10.1101/2022.12.21.521521v1)
  
  - DNA sequence embedding
    - [GROVER](https://huggingface.co/PoetschLab/GROVER)
    - [Model article](https://www.nature.com/articles/s42256-024-00872-0)
    model

- Manually trained model
  - Chaos Game Representation (CGR) and Frequency CGR (FCGR)
  - [Article](https://www.sciencedirect.com/science/article/pii/S2001037021004736)


