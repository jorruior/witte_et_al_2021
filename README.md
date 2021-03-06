# witte_et_al_2021
This folder contains the scripts used to generate the main output tables in the article 'A trans locus causes a ribosomopathy in hypertrophic hearts that affects mRNA translation in a protein length-dependent fashion' by Franziska Witte, Jorge Ruiz-Orera et al., currently accepted for publication in Genome Biology (2021). These scripts can be used to map RNA-seq and Ribo-seq datasets, and use these mapped files to predict QTLs and estimate heritability, as well as to call different categories of translated ORFs. With the generated output, we generated all Figures and Results. For the rest of methods, the used public software and parameters are described in the Methods section. License: GNU General Public License


**Mapping scripts:**

-*map_reads_TOPHAT.sh*: First bash script to map sequencing reads based on the Tophat mapper [https://doi.org/10.1186/gb-2013-14-4-r36]. Input parameters:
```
Argument 1: Fastq file(s) (single or paired). Filenames should end in '.fastq.gz'
Argument 2: GTF file. In our article, we used a modified rat GTF file based on the Ensembl v.82 release, including Ttn (not annotated) and masking the duplicated SURF locus (chr3:4,861,753-4,876,317 was masked and chr3:5,459,480-5,459,627 was included).
Argument 3: Tophat2 index to map the reads.
Argument 4: Bowtie2 index with a list of rRNA, tRNA and mtDNA contaminants that should be removed from the Ribo-seq datasets.
Argument 5: Read class. Arguments: 'single' (single-end RNA-seq), 'paired' (paired-end RNA-seq), 'ribo' (single-end Ribo-seq), or 'rna29' (single-end RNA-seq trimmed to 29bp).
Argument 6: Required for htseq-count. Read orientation: reverse or yes (forward).
```


-*map_reads_STAR.sh*: Second bash script to map sequencing reads based on the STAR mapper [https://doi.org/10.1093/bioinformatics/bts635]. Input parameters:
```
Argument 1: Fastq file(s) (single or paired). Filenames should end in '.fastq.gz'
Argument 2: GTF file. In our article, we used a modified rat GTF file based on the Ensembl v.82 release, including Ttn (not annotated) and masking the duplicated SURF locus (chr3:4,861,753-4,876,317 was masked and chr3:5,459,480-5,459,627 was included).
Argument 3: STAR index to map the reads.
Argument 4: Bowtie2 index with a list of rRNA, tRNA and mtDNA contaminants that should be removed from the Ribo-seq datasets.
Argument 5: Read class. Arguments: 'single' (single-end RNA-seq), 'paired' (paired-end RNA-seq), 'ribo' (single-end Ribo-seq), or 'rna29' (single-end RNA-seq trimmed to 29bp).
Argument 6: Required for htseq-count. Read orientation: reverse or yes (forward).
```


**Scripts for ORF prediction:**

To define the set of translated genes in rat heart and liver, we used RiboTaper v1.3 (Calviello et al., 2016) with standard settings  to  detect  open  reading  frames. Further processing of RiboTaper output using the scripts in *ORF_prediction*

-*select_ORFs_liver.R* and *select_ORFs_lv.R*: R script to extract different ORF types, e.g. CDS of known protein-coding genes, CDS of long non-coding RNAs, upstream ORFs etc. 

-*process_detected_ORFs.R*: R script to combine all detected ORFs and calculate FPKM in order to filter by minimum expression cutoff

**Scripts for QTL prediction:**

QTL mapping was performed using the linear regression model-based Matrix eQTL v2.1.1 (Shabalin, 2012) - scripts in *QTL_mapping*

-*1_matrixETL.R*: main R script to run matrixEQTL for the different datasets

-*2a_sort_qtl_results.R*, *2b_find_correct_SDP.R*, *2c_define_SDP_position.R* ad *2d_addFDR.R*: R scripts to further process QTL results

-*3_run_permutations.R*: R script to run permutation testing for all datasets

-*4_summarize_results_permutations.R*: R script to summarize all results and make venn diagrams

**Scripts to explore QTL results:**

-*tissue_comparison.R*


**Scripts for comparitive analysis of congenic rats:**

To  replicate the translatome-wide phenotype, we performed  ribosome  profiling  on  two  congenic  rat  lines (McDermott-Roe  et  al.,  2011). Scripts for comparative analysis are in *ComparativeAnalysis_CongenicRats* 

-*DESeq.R*

-*Assess_LengthEffect.R*

-*FC_GO_Term_plot.R*

**Other scripts:**

-*calculate_h2.R*: R script to estimate heritability based on the variance of replicates (BXH13 and BXH12). The script uses as input a table with the raw counts per replicate (3 vs 3, "replicate_counts.txt"). Afterwards, power was calculated using this online tool: https://github.com/Dashbrook/BXD_power_calculator_app

-*stoichometry.R*: R script to calculate and plot stoichometries based on the congenics and RI RNA-seq and Ribo-seq data. In this script, two genes are plotted: ENSRNOG00000008536 and ENSRNOG00000033734.

-*Process_CountData.R*: R script that contains a number of functions for NGS data wrangling, e.g. transformation of GTF into GenomicRanges object

**Required software**: 

-trim_galore (v.0.6.1): https://github.com/FelixKrueger/TrimGalore

-bowtie2 (v2.3.4.3): https://github.com/BenLangmead/bowtie2

-Tophat2 (v.2.1.0): https://ccb.jhu.edu/software/tophat/index.shtml

-STAR (v2.7.1a): https://github.com/alexdobin/STAR

-htseq-count (0.11.1): https://github.com/htseq/htseq
