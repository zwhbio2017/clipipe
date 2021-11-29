# CLIPipe

CLIPipe is an integrated pipeline for analyzing CLIP sequencing data. It provides all the commands needed to process CLIP-seq data, and it could identify and quantify sites of protein-RNA interactions on RNA from CLIP-seq data.

CLIPipe's input is the raw FASTQ file of CLIP-seq, and its main output is BED files with identified cross-linked sites.

Several other analyses are also included in CLIPipe that provide insights into the properties of protein-RNA interaction.

![Pipeline of Tutorial](img/CLIPipe_pipeline.png)

The CLIPipe workflow consists of:

-   Pre-processing function:
    -   Quality control, adapter removal, low-quality reads filtering, duplicates collapsing, and barcode removal of the raw CLIP-seq data.
-   Alignment function:
    -   Mapping pre-processed data to reference genome using bowtie, bwa, and novoalign
-   Peak calling function:
    -   Binding peak enrichment using Piranha, CTK, PureCLIP, iCLIPro, iCount, JAMM, PeakRanger, and clipcontext
-   Motif discovery function:
    -   Motif discovery of the binding regions using HOMER, PhyloGibbs, MEME, GraphProt, DREME and STREME

## Table of Contents:

-   [Requirements](1_requirement.md)
-   [Installation](2_installation.md)
-   [Basic Usage](3_basic_usage.md)
-   [Usage](4_usage.md)
    -   [Pre-processing](4_usage.md#pre-processing)
    -   [Alignment](4_usage.md#alignment)
    -   [Peak calling](4_usage.md#peak-calling)
    -   [Motif discovery](4_usage.md#motif-discovery)
-   [Frequently asked questions](5_faq.md)
-   [Citation](#citation)
-   [Copyright and License](6_copyright_and_license.md)

This documentation is for CLIPipe version: 1.0.3
