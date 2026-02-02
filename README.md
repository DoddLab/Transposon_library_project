# R Script for Figure Generation (Transposon library project)

## Overview
This repository contains **example R scripts** used for figure generation in the transposon library project. The scripts are provided for **methodological illustration** and can be adapted to new datasets with appropriate input files and paths. A detailed explanation of the scripts is provided below.



## Requirements
- R (>= 4.1 recommended)
- R packages used in the scripts (e.g., `tidyverse`, `ggplot2`, `ggtree`, `ggnewscale`, `VennDiagram`)

## Usage
1. Update file paths inside each script to point to your local data files.
2. Run the scripts in R or RStudio.
3. Outputs (figures and CSVs) are written to the locations specified in each script.


## Citation
- K. Sekiba, B. Hou, H. Chen, Z. Zhou, Y. Liu, A. Crawl-Bey, and D. Dodd*, High-throughput mass spectrometry reveals genetic determinants of nutrient acquisition in the gut bacterium Clostridium sporogenes, ***Nature Microbiology***, In Revision, 2026.


## Contact
For questions or issues, please open an issue on GitHub or contact Dylan Dodd or Kazuma Sekiba.


## License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a> 
This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)



---
# Script: `R_Figure_2D_Heat_map.R`

## Purpose
This file serves as a **test / example analysis script** demonstrating how second-round metabolite screening data from a transposon mutant library can be processed, analyzed, and visualized.

The script is provided for **figure generation and methodological illustration only**, and is not intended to represent a finalized or fully generalized analysis pipeline.


## Analysis Overview
The analysis integrates metabolite profiling data from multiple screening batches and links **metabolic phenotypes** to **genomic transposon insertion sites**.

The major goals are to:
- Quantify metabolite changes in individual mutants relative to WT
- Normalize effects using a biologically relevant reference strain (SACC)
- Identify mutants with statistically significant metabolic perturbations
- Visualize results in both metabolite-centric and genome-centric formats

---
## Input Data
The script uses multiple CSV files corresponding to independent screening batches (B1–B5).

Each input file contains:
- Mutant identifier (TLnumber)
- Gene annotation
- Transposon insertion position
- Quantified metabolite abundances

All datasets are reshaped into long format and combined prior to analysis.

---
## Key Analysis Steps
### 1. Data Preprocessing
- Convert wide metabolite tables to long format
- Replace missing values with zeros
- Exclude predefined metabolites not used for second screening
- Harmonize compound naming across batches

---
### 2. Mutant vs WT Comparison
For each metabolite:
- WT mean abundance is used as the baseline
- Each mutant is compared against WT to calculate:
  - log2 fold change
  - Absolute change
  - Relative change normalized by WT–SACC difference

Directionality of relative change is adjusted depending on metabolite class:
- **Substrates:** positive values indicate impaired utilization
- **Products:** positive values indicate reduced production

Statistical significance is assessed using two-sample t-tests where applicable.

---
### 3. Batch Integration and Filtering

- Results from all batches are merged
- Extreme values (Inf / NaN) are capped for visualization stability
- A curated list of valid mutants is applied
- Known low-quality or excluded mutants are removed

---
### 4. Heatmap-Style Phenotype Visualization
A dot-based heatmap summarizes results across mutants and metabolites:
- X-axis: transposon mutants (TLnumber, gene, insertion site)
- Y-axis: metabolites
- Color: normalized relative metabolic change
- Point size: statistical significance (–log10 p-value)

This figure enables rapid identification of:
- Pathway-specific metabolic defects
- Mutants with broad versus selective effects
- Consistency across metabolite classes

---
### 5. Genome-Level Visualization
Two genome-oriented figures are generated:
1. A genome-wide map of transposon insertion positions
2. A mapping of unique insertion sites to ordered phenotype “heat” coordinates

These figures facilitate genotype–phenotype interpretation at the genome scale.

---
## Output
- Filtered results table exported as CSV
- ggplot2-based figures suitable for publication or exploratory analysis

---
## Intended Use
This file is provided as a **test/example script** to:
- Explain analysis aims
- Illustrate data processing logic
- Demonstrate figure-generation strategies

Users are encouraged to adapt and extend this workflow for their own datasets.

<br><br>

---
# Script: `R_Figure_6M_Tree.R`

## Purpose
This file serves as a **test / example figure-generation script** demonstrating how to:

1) visualize a circular phylogenetic tree of DoddLab library strains,  
2) annotate the tree with **taxonomic metadata (Phylum)**,  
3) overlay multiple **per-strain continuous heatmaps** (e.g., BLASTP AA identity), and  
4) summarize gene/phenotype presence patterns using a **3-set Venn diagram**.

This example is intended for **methodological illustration and figure generation** rather than being a fully generalized pipeline.

## Overview of Analysis Aims
The script is designed to answer the following questions:

- How are strains distributed on a phylogeny, and what are their major taxonomic groupings?
- Which strains contain homologs of specific target genes (larA/B/C/E), and how strong is the similarity (AA identity)?
- How do these gene-level patterns overlap with a phenotype/marker readout (NPN.OCD)?
- Which strains satisfy combinations of these criteria?

## Inputs

### 1) Phylogenetic tree
- Newick tree file (protein-sequence consensus tree)
- Tip labels are cleaned to remove single quotes

### 2) BLASTP results table
- Long-format BLASTP output containing:
  - `strain`
  - `gene`
  - `AA.identity`
- The table is converted to wide format so each gene becomes a column.

### 3) Taxonomy metadata
- A mapping table linking strains (`tip.label`) to Phylum.

---
## Key Processing Steps
### 1. Tree tip label cleaning
- Removes `'` characters from `tree$tip.label` to standardize naming and enable merging.

### 2. BLASTP table reshaping and merge preparation
- Converts BLASTP output into wide format (`pivot_wider(names_from = gene, values_from = AA.identity)`).
- Renames `strain` → `tip.label` for compatibility with ggtree.
- Identifies strains present in the tree but missing from the BLASTP table, and appends them so the merge succeeds.

### 3. Missing-value handling and strain-level flags
- Replaces NA values with 0 for downstream filtering and plotting consistency.
- Defines a highlight label (`strain`) when `NPN.OCD >= 0.1`.
- Adds an `in_tree` flag and keeps only strains present in the phylogeny.

### 4. Label numbering for compact circular labeling
- Creates a numeric label (1…N) for each tree tip.
- Joins these labels to the metadata table used for plotting.
- Exports the final merged dataset to CSV for recordkeeping / manuscript use.

---
## Figure Generation
### Figure 1. Circular phylogeny annotated by Phylum
- Uses `ggtree(..., layout="circular")`.
- Tip points are colored by **Phylum**.
- Tip labels use compact numeric identifiers (`label_number`) to reduce clutter.
- Strains meeting a threshold (e.g., `NPN.OCD >= 0.1`) are additionally labeled in red.

### Figure 2. Multi-track heatmaps aligned to the phylogeny
Heatmaps are added sequentially using `gheatmap()` + `ggnewscale::new_scale_fill()`:

- larB AA identity
- larC AA identity
- larE AA identity
- larA AA identity
- NPN.OCD AA identity

Each heatmap track:
- is aligned to the same tree order,
- uses a continuous fill scale,
- uses defined breaks and limits for interpretability,
- sets `na.value = "white"` to cleanly display missing data.

---
## Venn Diagram Summary

A 3-set Venn diagram summarizes overlaps across strain sets defined by thresholds:
- Group1: `larB >= 0.3 & larC >= 0.3 & larE >= 0.3`
- Group2: `larA >= 0.3`
- Group3: `NPN.OCD >= 0.3`

This visualization provides a compact overview of:
- shared strain membership,
- unique strain subsets,
- candidate strains meeting multiple criteria simultaneously.

---
## Outputs
- `For_paper.csv`: merged per-strain table including
  - AA identity per gene
  - NPN.OCD values
  - taxonomy (Phylum)
  - tip label numbering
- ggplot2-based figures:
  - circular phylogeny with labels
  - heatmap tracks aligned to the tree
- Venn diagram (3-set overlap)

---
## Notes / Intended Use

- File paths are currently set to local directories and should be updated for reproducible demos.
- Thresholds (0.1, 0.3) are example cutoffs chosen for exploratory visualization and can be tuned.
- This file is meant as a **test/example** for figure generation and downstream manuscript preparation.

---