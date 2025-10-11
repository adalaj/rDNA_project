
## In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus(2025) script details:

This folder contains all the scripts for identification of predicted non-canonical DNA structure sequences in chicken rDNA sequence.



## Requirements

The script requires the following R libraries:
 - `data.table`
 - `tidyverse`
 - `karyoplote`
 - `Biostrings`
 

**Note:** Ensure the required packages are installed before running the script 


## Chicken Folder

The upstream pipelines used to generate these files (sequence scanning and motif prediction) are described in the Methods section of the paper and are not reproduced here. The files are included here only as inputs for figure generation.


## Scripts and descriptions

#### 1. `chicken_evolutionary.R`
**Inputs:**

- `G4FS_KT445934_chicken_junctn_details.csv`

- `RLFS_KT445934_chicken_junctn_details.csv`

- `chicken_imotif_prediction_end_to_end_prediction_default_setting_master.csv`

These files contain genomic coordinates for G4-forming sequences (G4FS), R-loop forming sequences (RLFS), and i-motif forming sequences (iMFS). All are deposited in supplementary tables.

**Purpose:** Integrates these datasets for comparative and evolutionary visualization.

---

#### 2. `chicken_imotif.R`
**Input:** iM-seeker (https://im-seeker.org/) output CSV files for template and non-template strands  
**Output:** `chicken_imotif_prediction_end_to_end_prediction_default_setting_master.csv`  
**Purpose:** Assigns iMFS predictions to defined rDNA subregions and generates summary outputs.

---

#### 3. `chicken_rdna_details.R`
**Input:** Chicken rDNA FASTA sequence (GenBank: *KT445934*)(https://www.ncbi.nlm.nih.gov/nuccore/KT445934)
**Output:** Nucleotide composition and GC-content distribution CSV file  
**Purpose:** Computes base composition metrics across the rDNA locus.

---

#### 4. `chicken_rdna_G4FS.R`
**Input:** Output BED file from `g4_canonical_finder_3.11python.py` (https://github.com/CsabaPapp13/Stable-bulged-G-quadruplexes-in-the-human-genome_2022-2023)
**Output:** `G4FS_KT445934_chicken_junctn_details.csv`  
**Purpose:** Assigns G4-forming sequences to rDNA subregions, normalizes counts, and produces distribution plots.

---

#### 5. `chicken_rdna_rloop.R`
**Input:** Output BED file from QmRLFS-finder (v1.5) (http://r-loop.org/?pg=qmrlfs)  
**Output:** `RLFS_KT445934_chicken_junctn_details.csv`  
**Purpose:** Maps predicted R-loop forming sequences to rDNA components, visualizes strand distributions, and generates normalized plots.

---


## Contact

For any questions or issues, contact Jyoti Devendra Adala.
