# 🧬 E. coli Transcriptional Regulatory Network Analyzer

An advanced, modular Python toolkit for analyzing the Transcriptional Regulatory Network (TRN) of *Escherichia coli*. 

This project demonstrates the application of **Network Science** and **Bioinformatics** to model complex biological systems. It extracts genomic features from GenBank files, infers operon structures based on intergenic distances, and builds a directed graph of Transcription Factors (TFs) and their Target Genes (TGs).

## ✨ Features & Capabilities

* **Genomic Parsing (`parsing_utils.py`):** Utilizes `Biopython` to parse GenBank annotations (`.gb`) and map TSV regulatory interactions into a robust data structure. Includes an algorithm for operon prediction based on strand orientation and spatial constraints.
* **Topological Analysis (`network_metrics.py`):** Calculates key network metrics using `networkx`, including degree distributions (in-degree/out-degree hubs), Largest Connected Components (LCC), and graph diameter estimations via randomized pathfinding.
* **Network Motifs & Statistical Significance (`motifs_analysis.py`):** Implements algorithms to discover critical structural motifs like Feed-Forward Loops (FFL). It features a robust validation engine that generates null models (Randomized Configuration Models) to statistically prove the evolutionary significance of these motifs against background noise.
* **Modular Pipeline (`main.py`):** A clean, execution-ready orchestrator script that ties all modules together.

## 🛠️ Tech Stack & Requirements

* **Python 3.8+**
* `NetworkX`: Graph creation and topology algorithms.
* `Biopython`: Parsing biological sequence records.
* `NumPy`: Statistical calculations (mean, standard deviations).

Install the requirements via pip:
```bash
pip install networkx biopython numpy
```

## 🚀 How to Run

1. Clone this repository.
2. Ensure you have the required data files in the root directory:
   * `U00096.3.gb` (E. coli K-12 complete genome)
   * `TFSet.tsv` (Transcription factors dataset)
   * `TF-RISet.tsv` (Regulatory interactions dataset)
3. Execute the main pipeline:
   ```bash
   python main.py
   ```
4. The script will output topological statistics to the terminal and generate two `.graphml` files (`Ecoli_TRN.graphml` and `Ecoli_operon_TRN.graphml`) that can be visualized in software like **Gephi** or **Cytoscape**.

## 🏗️ Project Architecture

The codebase strictly follows the *Separation of Concerns* principle, making it scalable for future bioinformatics pipelines. Data extraction is decoupled from mathematical graph analysis, which is in turn decoupled from statistical modeling.
