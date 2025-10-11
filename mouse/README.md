
## In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus(2025) script details:


##Mouse Folder
This folder contains all the scripts for identification of predicted non-canonical DNA structure sequences in mouse rDNA sequence.



## Requirements

The script requires the following R libraries:
 - `data.table`
 - `tidyverse`
 - `karyoplote`
 - `Biostrings`
 - `seqinr`
 - `stringr`
 

**Note:** Ensure the required packages are installed before running the script 


## Scripts and descriptions

#### 1. `mouse_evolutionary.R`
**Inputs:**

- `G4FS_BK000964_added_5000nt_IGS_upstream_at_junctn_details.csv`

- `RLFS_BK000964_added_5000nt_IGS_upstream_master_qmrlfs_table_after_rule.csv`

- `mouse_imotif_prediction_end_to_end_prediction_default_setting_master.csv`

These files contain genomic coordinates for G4-forming sequences (G4FS), R-loop forming sequences (RLFS), and i-motif forming sequences (iMFS). All are deposited in supplementary tables 8, 9 and 10 respectively.

**Purpose:** Integrates these datasets for comparative and evolutionary visualization.

---

#### 2. `mouse_imotif.R`
**Input:** iM-seeker (https://im-seeker.org/) output CSV files for template and non-template strands  
**Output:** `mouse_imotif_prediction_end_to_end_prediction_default_setting_master.csv`  
**Purpose:** Assigns iMFS predictions to defined rDNA subregions and generates summary outputs.

---

#### 3. `mouse_rdna_details.R`
**Input:** mouse rDNA FASTA sequence (GenBank: *BK000964*)(https://www.ncbi.nlm.nih.gov/nuccore/BK000964?report=genbank)

**Output:** Nucleotide composition and GC-content distribution CSV file  
**Purpose:** Computes base composition metrics across the rDNA locus.

---

#### 4. `mouse_rdna_G4FS.R`
**Input:** Output BED file from `g4_canonical_finder_3.11python.py` (https://github.com/CsabaPapp13/Stable-bulged-G-quadruplexes-in-the-human-genome_2022-2023)

**Output:** `G4FS_BK000964_added_5000nt_IGS_upstream_at_junctn_details.csv`  
**Purpose:** Assigns G4-forming sequences to rDNA subregions, normalizes counts, and produces distribution plots.

---

#### 5. `mouse_rdna_rloop.R`
**Input:** Output BED file from QmRLFS-finder (v1.5) (http://r-loop.org/?pg=qmrlfs)  
**Output:** `RLFS_BK000964_added_5000nt_IGS_upstream_master_qmrlfs_table_after_rule.csv`  
**Purpose:** Maps predicted R-loop forming sequences to rDNA components, visualizes strand distributions, and generates normalized plots.

---


## Contact

For any questions or issues, contact Jyoti Devendra Adala.
