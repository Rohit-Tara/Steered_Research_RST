# Steered_Research_RST 
# Replication of scRNA-seq Pipeline  
## Replication of Camp et al. (2015) Figure 3D  
### University of Leicester — BS7120 Steered Research Project Coursework  

---

## Overview

This repository contains the **original RNA-seq pipeline** used to replicate Figure 3D from:

> Camp, J.G. et al. (2015). *Human cerebral organoids recapitulate gene expression programs of fetal neocortex development*.  
> PNAS, 112(51), 15672–15677.  
> https://doi.org/10.1073/pnas.1520760112  

The original paper used **TopHat2 for alignment** and **Cufflinks for transcript quantification**. This pipeline closely follows that approach, with downstream clustering and visualisation performed using **Seurat**.

---

## Data

Raw sequencing data: GEO accession [GSE75140](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75140)

- 734 single cells (508 organoid, 226 fetal neocortex)  
- Paired-end reads   
- Illumina   

Reference genome: **GENCODE GRCh38 release 22**

- Genome FASTA: `GRCh38.primary_assembly.genome.fa`  
- GTF annotation: `gencode.v22.primary_assembly.annotation.gtf`  
- https://www.gencodegenes.org/human/release_22.html  

---

## Scripts 
Run in order of top to bottom:
| Script | Description |
|--------|------------|
| `fastq_convert.sh` | Convert the downloaded SRA files to FASTQ format |
| `retry_missing_fastq.sh` | Retry failed or incomplete FASTQ conversions |
| `nochr.py` | Remove chromosome prefixes from GTF for compatibility |
| `tophat_all2.slurm` | Align reads to GRCh38 using TopHat2 (HPC job script) |
| `cufflinks_top.sh` | Quantify gene expression (FPKM values) |
| `expression_matrix.py` | Merge Cufflinks outputs into a unified expression matrix |
| `v5_steered.R` | Seurat analysis (PCA, clustering, t-SNE, marker genes) | 

## Requirements

- HPC cluster with SLURM job scheduler for alignment step
- Minimum 32GB RAM recommended for Seurat analysis
- R packages: Seurat, ggplot2, dplyr, data.table, ggdendro, patchwork

Install R packages with:
```r
install.packages(c("Seurat", "ggplot2", "dplyr", 
                   "data.table", "ggdendro", "patchwork"))
```

---

## Software Versions
| Tool | Version |
|------|---------|
| SRA Toolkit | 3.0.7 |
| TopHat2 | 2.1.1 |
| Cufflinks | 2.2.1 |
| Bowtie2 | 2.4.4 |
| Python | 3.9 |
| R | 4.5.2 |
| Seurat | 5.0 | 

## Usage

1. Download raw data and the metadata `SraRunTable.csv` from SRA using accession SRP066834:(https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA304502&o=acc_s%3Aa)
2. Run scripts in order as listed in the Scripts table above
3. Ensure reference genome and GTF are downloaded before running `tophat_all2.slurm`
4. Output expression matrix will be saved as `expression_matrix_full.csv`
5. Load `expression_matrix_full.csv` and `SraRunTable.csv` into `v5_steered.R` for clustering analysis

## Citation 
Camp, J.G. et al. (2015) 'Human cerebral organoids recapitulate gene expression programs of fetal neocortex development', Proceedings of the National Academy of Sciences of the United States of America, 112(51), pp. 15672–15677. Available at: https://doi.org/10.1073/pnas.1520760112 . 
