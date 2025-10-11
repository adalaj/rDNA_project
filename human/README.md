
## In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus(2025) script details:

## Human Folder
This folder contains all the scripts and graph input files needed for identification of predicted non-canonical DNA structure sequences in human rDNA sequence. Each script descriptions are listed below. The archive folder contains all the scripts that are either not needed or are incomplete.  


## Requirements

The script requires the following R libraries:
 - `data.table`
 - `tidyverse`
 - `karyoploteR`
 - `Biostrings`
 - `corrplot`
 - `VennDiagram`
 

**Note:** 
 - Ensure the required packages are installed before running the script.
 - All output files are provided either in the supplementary tables or within this directory.

---

## Scripts and descriptions


#### 1. `G4FS_at_human_rDNA.R`
**Input:** Output BED file from `g4_canonical_finder_3.11python.py` (https://github.com/CsabaPapp13/Stable-bulged-G-quadruplexes-in-the-human-genome_2022-2023)

**Output:** `G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv`,  deposited in supplementary table 3. 
This output file was subsequently processed to generate the input and output files for figure 3C, 3D and 3E. 
Figure 3C and 3D use `G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv`, while Fig 3E use `G4FS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_graphinput.csv` graph input files.

**Purpose:** Assigns G4-forming sequences to rDNA subregions, normalizes counts, and produces distribution plots.

---

#### 2. `G4FS_rdna_locus_visualization.R`
**Input:** Output BED file from `g4_canonical_finder_3.11python.py` (https://github.com/CsabaPapp13/Stable-bulged-G-quadruplexes-in-the-human-genome_2022-2023)

**Output:** Generated fig 3F and 3G for human.

**Purpose:** Visualisation of G4-forming sequences at human rDNA.

---


#### 3. `GC_content_function.R`
**Input:** Any DNA sequence.

**Output:** Generate GC content for the sequence. 

**Purpose:** To calculate GC content for any DNA sequence. 

---


#### 4. `gc_skew_function.R`
**Input:** Any DNA sequence.

**Output:** Generates GC skew values for the input DNA sequence.

**Purpose:** Calculates GC skew to assess strand compositional asymmetry.

---


#### 5. `human_evolutionary.R`
**Inputs:**

- `G4FS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv`

- `RLFS_KY962518_added_3500nt_IGS_upstream_master_qmrlfs_table_after_rule.csv`

- `imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv`

These files contain genomic coordinates for G4-forming sequences (G4FS), R-loop forming sequences (RLFS), and i-motif forming sequences (iMFS) in human rDNA sequence. 
All are deposited in supplementary tables 2, 3 and 4 respectively. These input files were filtered for sequences found only in human rDNA transcriptional unit (5'ETS to 3'ETS). 

**Output:** Generated fig 7 for Human.

**Purpose:** Integrates these datasets, filter them for human rDNA transcriptional unit (5'ETS to 3'ETS) for comparative and evolutionary visualization.

---


#### 6. `human_rdna_details.R`
**Input:** human rDNA FASTA sequence (GenBank: *KY962518*)(https://www.ncbi.nlm.nih.gov/nuccore/KY962518.1)

**Output:** CSV tables and bar plots showing nucleotide (A, T, G, C) and AT/GC content across defined rDNA subregions. Output file is depoisted in supplementary table 1.

**Purpose:** Extracts annotated rDNA domains, calculates nucleotide composition, and generates figure-ready visualizations (Fig 1B–D)

---


#### 7. `human_rdna_gc_skew.R`
**Input:** human rDNA FASTA sequence (GenBank: *KY962518*)(https://www.ncbi.nlm.nih.gov/nuccore/KY962518.1), `GC_content_function.R` and `gc_skew_function.R`

**Output:** Sliding-window CSV tables and plots of GC content (%) and GC skew across the rDNA locus.

**Purpose:** Quantifies and visualizes strand bias and base composition patterns (Fig 1E–F)

---


#### 8. `human_vs_chicken_rdna_global_similarity.R`
**Input:** Human (*KY962518*) and chicken (*KT445934*) rDNA FASTA files with region annotations.  
**Output:** CSV tables summarizing pairwise sequence identity across rDNA regions. Deposited in supplementary table 7.  
**Purpose:** Performs global alignments to quantify sequence conservation between human and chicken rDNA loci.

---


#### 9. `human_vs_monkey_rdna_similarity.R`
**Input:** Human (*KY962518*) and monkey (*KX061890*) rDNA FASTA files with region annotations.  
**Output:** CSV tables summarizing pairwise sequence identity across rDNA regions. Deposited in supplementary table 7.  
**Purpose:** Performs global alignments to quantify sequence conservation between human and monkey rDNA loci.

---

#### 10. `human_vs_mouse_rdna_similarity.R`
**Input:** Human (*KY962518*) and mouse (*BK000964*) rDNA FASTA files with region annotations.  
**Output:** CSV tables summarizing pairwise sequence identity across rDNA regions. Deposited in supplementary table 7.  
**Purpose:** Performs global alignments to quantify sequence conservation between human and mouse rDNA loci.

---


#### 11. `iMFS_at_human_rDNA.R`
**Input:** iM-seeker (https://im-seeker.org/) output CSV files for template and non-template strands  
**Output:** `imotif_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv`,  deposited in supplementary table 4. 
This output file was subsequently processed to generate the input and output files for figure 4C, 4D and 4E. 
Figure 4C and 4D use `iMFS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv`, while Fig 4E use `imotif_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_AR_graphinput.csv` graph input files.
Generated fig 4F and 4G for human.

**Purpose:** Assigns iM-forming sequences to rDNA subregions, normalizes counts,produces distribution plots, and Visualisation of iM-forming sequences at human rDNA.


---

#### 12. `noncanonical_vs_polr1_local_variation.R`
**Input:** ChIP-seq counts for POLR1A  with RLFS, G4FS, and iMFS annotations across *chrR* (hg38-rDNA_v1.0). 
ChIP seq and chrR was taken from (https://github.com/vikramparalkar/rDNA-Mapping-Genomes) git repository.

**Output:** Correlation plots, line graphs, and summary CSVs linking Pol I occupancy to non-canonical DNA structures.
Fig 6 and supplementary table 6.

**Purpose:** Assesses co- or anti-correlation between transcriptional machinery and predicted R-loops, G4s, and i-motifs along the rDNA locus.

---

#### 13. `RIZ_G4s_imotif_correlation.R`
**Input:** Processed coordinate files for RIZ, G4FS, and iMFS across human rDNA  
**Output:** CSV and plots showing binned distributions, correlation matrix, Venn overlap, and scatter relationships.
Refer fig5 and supplementary Table 5.

**Purpose:** Compares positional enrichment and co-occurrence of RIZ, G4FS, and iMFS within the rDNA locus.

---



#### 14. `RIZ_only_at_human_rdna.R`
**Input:** Output BED file from QmRLFS_RIZ_finder.py (modified from http://r-loop.org/?pg=qmrlfs)  

**Output:** `RIZ_KY962518_added_3500nt_IGS_upstream_at_junctn_details_after_rule.csv` deposited in repository. 

**Purpose:** Mapping of RIZ in human rDNA subregions

---

#### 15. `RLFS_at_human_rDNA.R`
**Input:** Output BED file from QmRLFS-finder (v1.5) (http://r-loop.org/?pg=qmrlfs)  

**Output:** `RLFS_KY962518_added_3500nt_IGS_upstream_at_junctn_details.csv`,  deposited in supplementary table 2. 
This output file was subsequently processed to generate the input and output files for figure 2C, 2D and 2E. 
Figure 2C and 2D use `RLFS_KY962518_added_3500nt_IGS_upstream_at_junctn_after_rule_graphinput.csv`, while Fig 2E use `RLFS_KY962518_added_3500nt_IGS_upstream_no_junctn_strandwise_graphinput.csv` graph input files.

**Purpose:** Assigns Rloop-forming sequences to rDNA subregions, normalizes counts, and produces distribution plots.

---


#### 16. `RLFS_rdna_locus_visualisation.R`
**Input:** Output BED file from QmRLFS-finder (v1.5) (http://r-loop.org/?pg=qmrlfs)  
**Output:** Generated fig 2F and 2G for human.
**Purpose:** Visualisation of Rloop-forming sequences at human rDNA.

---

#### 17. `QmRLFS_RIZ_finder.py`
**Input:** Any DNA sequence in fasta format.
**Output:** Provide information of genomic coordinates of R-loop initiation zone (RIZ) in strand specific manner.
**Purpose:** To find RIZ in human rDNA.

---


## Contact

For any questions or issues, contact Jyoti Devendra Adala.
