## In silico Mapping of Non-Canonical DNA Structures Across the Human Ribosomal DNA Locus(2025):
This repository contains all scripts used for in silico prediction, analysis, and visualization of non-canonical DNA structures (NCSs) — including R-loops (RLFS/RIZ), G-quadruplexes (G4FS), and i-motifs (iMFS) — across the ribosomal DNA (rDNA) loci of human, monkey, mouse, and chicken.
Each organism-specific folder contains the R scripts, auxiliary Python scripts, and graph input files required for figure generation and analysis described in the associated manuscript.


## Folder Structure
- `chicken`– Scripts and data for Gallus gallus rDNA locus
- `human` – Scripts and data for Homo sapiens rDNA locus
- `monkey`– Scripts and data for Macaca mulatta rDNA locus
- `mouse` – Scripts and data for Mus musculus rDNA locus

An `archive` subfolder contains incomplete or redundant scripts
 
Each species folder contains:
1. R scripts for GC-based metrics, RLFS, G4FS, and iMFS mapping
2. Evolutionary comparison scripts
3. Figure-ready CSV and visualization outputs
4. A detailed `README.md` files describing script puposes and input/outputs 


## Requirements

The script requires the following R libraries:
 - `data.table`
 - `tidyverse`
 - `karyoploteR`
 - `Biostrings`
 - `corrplot`
 - `VennDiagram`
 

**Note:** 
1. Ensure that all required packages are installed before running the scripts.
2. All output files are provided either in the supplementary tables or within their respective species folders.


## Contact

For any questions or issues, contact Jyoti Devendra Adala (adalaj@upstate.edu).
