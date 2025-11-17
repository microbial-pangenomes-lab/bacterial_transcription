# Cis non-coding genetic variation drives gene expression changes in the *E. coli* and *P. aeruginosa* pangenomes

## Overview of Analytical Notebooks & Scripts
Below is a short description for the scripts and notebooks used for all analyses presented in the manuscript. The descriptions below provide readers with a clear overview of the purpose and function of each script.

## Differential expression analyses

## *differential_analysis_DESeq2.R* ##
Performs differential expression analysis using DESeq2 on raw RNA-seq count matrices for *E. coli* and *P. aeruginosa* strains.

*differential_analysis_LPEseq.R*
Performs differential expression analysis using LPEseq on raw RNA-seq count matrices for *E. coli* and *P. aeruginosa* strains.

*DEGs_ecoli_amr.ipynb*
Conducts differential expression analysis for *E. coli* and *P. aeruginosa* strains with and without resistance-associated variants exposed to Tobramycin.

## Cis-regulatory variants and phenotype associations

*gene_expression_phenotype.ipynb*
Selects genes to be tested for associations with cis-regulatory variants

*cis-GWAS_gene_expression.ipynb* 
Analyses pyseer GWAS outputs to identify cis-regulatory variants associated with gene-expression variation.

*cis-GWAS_resistance_phenotype.ipynb*
Analyses pyseer GWAS outputs to identify cis-regulatory variants associated with antimicrobial resistance phenotypes.

*wholegenome-GWAS_AMR.ipynb*
Performs whole-genome association analyses using core variants, accessory genes, and k-mers to detect genome-wide determinants of antimicrobial resistance.

*validating_correlations.ipynb*
Validates variant–expression and expression–phenotype associations using independent statistical approaches and alternative correlation models.

*manhattan-qq_plots.ipynb*
Generates Manhattan and QQ plots from GWAS outputs to visualise genome-wide association signals and p-value distributions.

## Transcriptional Modeling and Regulatory Interpretation

*predicted_transcription_rates.ipynb*
Computes predicted transcription rates for each gene or operon using sequence-based features (promoter content, motifs, UTR features, etc.).

*nCV_TFBS_overlap.ipynb*
Assesses overlap between significant non-coding variants and transcription factor binding sites (TFBS) to infer regulatory mechanisms underlying cis-effects.

*aln2img_modified2.py*
Converts aligned nucleotide sequences into heatmap-style images, enabling visualisation of promoter diversity, motif patterns, and non-coding variant structure across strains.

## Experimental and Utility Scripts

*multi_day_gfp_analysis.ipynb*
Processes multi-day GFP reporter assay data to quantify promoter activity and compare regulatory constructs over time.

*chone's_test_foldenrichment.ipynb*
Computes fold-enrichment of variants or genes within selected genomic features or operons, providing support for enrichment-based analyses.
