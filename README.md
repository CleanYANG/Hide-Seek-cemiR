# Hide-Seek-cemiR

A small tool to predict **miRNA–RNA interactions** using:
- canonical seed types (8mer, 7mer-m8, 7mer-A1)
- RNA secondary structure (RNAfold) to estimate accessibility
- local duplex free energy (ΔG) using RNAduplex
- shared miRNAs between two RNAs (e.g. lncRNA Stgart and 3'UTR-Star)
- optional PubMed annotation

## Installation

```bash
mamba env create -f environment.yml
mamba activate hide-seek-cemir
