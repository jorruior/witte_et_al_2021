# witte_et_al_2021
This folder contains the scripts used to generate the main output tables in the article 'A trans locus causes a ribosomopathy in hypertrophic hearts that affects mRNA translation in a protein length-dependent fashion' by Franziska Witte, Jorge Ruiz-Orera et al., currently accepted for publication in Genome Biology (2021). These scripts can be used to map RNA-seq and Ribo-seq datasets, and use these mapped files to predict QTLs and estimate heritability, as well as to call different categories of translated ORFs. With the generated output, we generated all Figures and Results. For the rest of methods, the used public software and parameters are described in the Methods section. License: GNU General Public License


**Mapping scripts:**

--


-*map_reads_v2.sh*: Second bash script to map sequencing reads based on the STAR mapper [https://www.ncbi.nlm.nih.gov/pubmed/23104886]. Input parameters:
```
Argument 1: Fastq file(s) (single or paired). Filenames should end in '.fastq.gz'
Argument 2: GTF file. In our article, we used a modified rat GTF file based on the Ensembl v.82 release, including *Ttn* gene and masking the duplicated SURF locus (chr3:4,861,753-4,876,317 was masked and chr3:5,459,480-5,459,627 was included).
Argument 3: STAR index to map the reads.
Argument 4: Bowtie2 index with a list of rRNA, tRNA and mtDNA contaminants that should be removed from the Ribo-seq datasets.
Argument 5: Class of reads. single (single-end RNA-seq), paired (paired-end RNA-seq), or ribo (single-end Ribo-seq).
Argument 5: Required for htseq-count. Read orientation: reverse or yes (forward).
```


**Scripts for ORF prediction:**

To define the set of translated genes in rat heart and liver, we used RiboTaper v1.3 (Calviello et al., 2016) with standard settings  to  detect  open  reading  frames.

Further processing of RiboTaper output using the scripts in *ORF_Detection*

-*select_ORFs_liver.R* and *select_ORFs_lv.R*: R script to extract different ORF types, e.g. CDS of known protein-coding genes, CDS of long non-coding RNAs, upstream ORFs etc. 

-*process_detected_ORFs.R*: R script to combine all detected ORFs and calculate FPKM in order to filter by minimum expression cutoff

**Scripts for QTL prediction:**

--

**Other scripts:**

-*calculate_h2.R*: R script to estimate heritability based on the variance of replicates (BXH13 and BXH12). The script uses as input a table with the raw counts per replicate (3 vs 3, "replicate_counts.txt"). Afterwards, power was calculated using this online tool: https://github.com/Dashbrook/BXD_power_calculator_app

-*stoichometry.R*: R script to calculate and plot stoichometries based on the congenics and RI RNA-seq and Ribo-seq data. In this script, two genes are plotted: ENSRNOG00000008536 and ENSRNOG00000033734.



**Required software**: 

-trim_galore (v.0.6.1): https://github.com/FelixKrueger/TrimGalore

-bowtie2 (v2.3.4.3): https://github.com/BenLangmead/bowtie2

-STAR (v2.7.1a): https://github.com/alexdobin/STAR

-htseq-count (0.11.1): https://github.com/htseq/htseq
