
# CLIPipe

CLIPipe is an integrated pipeline for analyzing CLIP sequencing data. It provides all the commands needed to process CLIP-seq data, and it could identify and quantify sites of protein-RNA interactions on RNA from CLIP-seq data.

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

- [Requirements](#requirements)
- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Usage](#usage)
  - [Pre-processing](#pre-processing)
  - [Alignment](#alignment)
  - [Peak calling](#peak-calling)
  - [Motif discovery](#motif-discovery)
- [Copyright and License Information](#copyright-and-license-information)
- [Citation](#citation)
- [Tutorial](https://clipipe.readthedocs.io/)

## Requirements

### Software
All required software and packages are already installed in the docker, so there are no more requirements. You just need to install docker.

-   docker


## Installation

### Use CLIPipe Docker

#### Prepare reference demo and docker

1. Download the reference, demo data, and CLIPipe docker to directory, and unzip

         mkdir clipipe_test;
         cd clipipe_test;
         wget http://clipipe.ncrnalab.org/clipipe_ref.tar.gz;
         wget http://clipipe.ncrnalab.org/clipipe_demo.tar.gz;
         wget http://clipipe.ncrnalab.org/CLIPipe_v1.0.3_.tar.gz;
         tar -xvzf clipipe_ref.tar.gz;
         tar -xvzf clipipe_demo.tar.gz;

#### Run CLIPipe Docker

1.  [Docker](https://www.docker.com/) provides an easy way to run CLIPipe in a working environment that is completely separated from your host machine. All required software and packages are already installed in a ready-to-use image of CLIPipe docker, so there are no more requirements. You can use the docker image we provide: [CLIPipe Docker Image](https://hub.docker.com/repository/docker/shangzhang/clipipe). And you may execute these commands to get the docker `CLIPipe_1.0.X` container:
         
         cd clipipe_test;
         docker import CLIPipe_v1.0.3_.tar.gz zs/clipipe:1.0.3_test     ##import the docker

         docker run --name=CLIPipe_1.0.3_test -t -d -h CLIPipe_docker --restart unless-stopped -v <the-absolute-path-of-current-directory>:/home/CLIPipe_user/clipipe zs/clipipe:1.0.3_test /bin/bash

    -   Make sure to create a local folder and provide the path to it. The example above uses a path that may not be applicable to your setup. Both, path to the folder on the host machine and path within the container (`/home/CLIPipe_user/clipipe`) must be absolute.

2.  To show the docker container `CLIPipe_1.0.3_test`, you can execute:

         docker container ls

3.  To execute the contianer `CLIPipe_1.0.3_test`, you can execute:

         docker exec -it CLIPipe_1.0.3_test bash

         cd /home/CLIPipe_user/clipipe;
         find clipipe_demo  -type d  -exec chmod 777 {} \;
         find clipipe_demo  -type f -exec chmod 666 {} \;

4.  After entering the container, please change the user to `CLIPipe_user`

         su CLIPipe_user;
         cd ~

5.  To test the installation and get information about the command-line interface of CLIPipe, you can execute:

         clipipe --help

    A helper message is shown like this:

         usage: clipipe [-h] --user_config_file USER_CONFIG_FILE
                        {pre_process,mapping,peak_calling}

         CLIPipe: A comprehensive quality control and analysis pipeline for CLIP highthroughput sequencing data
         =======================================================================================================
         CLIPipe is a Python module and associated command-line interface (CLI), which provides all the
         commands needed to process protein-RNA CLIP interaction data and to identify and quantify
         sites of protein-RNA interactions on RNA.

         CLIPipe's main input are FASTQ files with CLIP-seq data, its main output are BED files
         with identified and quantified cross-linked sites.

         A number of analyses are included in CLIPipe that provide insights into the properties of
         protein-RNA interaction.

         optional arguments:
            -h, --help                show this help message and exit
            -u, --user_config_file    the user config file

         positional arguments:
            {pre_process,mapping,peak_calling}

         =======================================================================================================
         For additional help or support, please visit https://github.com/ShangZhang/clipipe

### Use from source

We have hide the details of each step using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and you only need to run one single command. However, you can use some of the codes if you are familiar with snakemake. The source code is [here](https://github.com/ShangZhang/clipipe).

## Basic Usage
The basic usage of CLIPipe is:

```bash
clipipe ${step_name} -d ${dataset}

# Note:
#   ${step_name} is one of the step listed in 'positional arguments'.
#   ${dataset} is the name of your dataset that should match the prefix of your configuration file described in the following section.
```

## Usage

### Reference data

You could use the provided reference file to run CLIPipe. Defaultly you may choose from hg19/hg38 (human) and mm10/mm39 (mouse). You can also create your own reference based on the scripts that we provided.

```bash
ls /home/CLIPipe_user/clipipe/clipipe_ref/
```

### Demo data

You can use the provided demo data to run CLIPipe:

```bash
cd /home/CLIPipe_user/clipipe/clipipe_demo/general/
```

The demo data folder has the following structure:

```text
./
├── config
|   ├── default_config.yaml
│   └── user_config.yaml
├── data
|   ├── fastq/
│   └── sample_ids.txt
└── output
    └── ...
```

```text
Note:
    `config/user_config.yaml`: configuration file with user defined parameters for each step.
    `config/default_config.yaml`: configuration file with additional detailed parameters for each step. The default file is not supposed to be changed unless you are very clear about what you are doing.
    `config/fastq/`: folder of raw CLIP-seq fastq file.
    `data/sample_ids.txt`: table of sample name information.
    `output/example/`: output folder.
```

### User config file
The user config file is shown like this:
```text
# default config
default_config_file: /home/CLIPipe_user/clipipe/clipipe_demo/general/config/default_config.yaml

# basic config file path
species: Human_hg38
reference_dir: /home/CLIPipe_user/clipipe/clipipe_ref

data_dir: /home/CLIPipe_user/clipipe/clipipe_demo/general/data
temp_dir: /home/CLIPipe_user/clipipe/clipipe_demo/general/temp
output_dir: /home/CLIPipe_user/clipipe/clipipe_demo/general/output_human_hg38
summary_dir: /home/CLIPipe_user/clipipe/clipipe_demo/general/summary

# general parameters
threads_compress: 2
threads_mapping: 4

# pre process parameters
paired_end: false
barcode_length: 1
adaptor1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
adaptor2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# mapping parameters
aligner: bowtie

# peak calling parameters
peak_caller: Piranha
read_length: 100
```

### Pre-processing

CLIPipe provides pre-process step for raw CLIP-seq data. You need to set up the `config/user_config.yaml` file correctly. The other parameters for pre-process step can be found in `config/default_config.yaml`.

```bash
cd /home/CLIPipe_user/clipipe/clipipe_demo/general/;
clipipe -u ./config/user_config.yaml pre_process
```

```text
Note:
    The output folder `output_human_hg38/fastqc_raw/` contains quality control results of raw CLIP-seq data.
    The output folder `output_human_hg38/multiqc_raw/` contains summary of all raw sequencing data quality control results.
    The output folders `output_human_hg38/pre_process/` contain the pre process results of raw CLIP-seq data.
```

### Alignment

CLIPipe provides bowtie, bwa and novoalign for mapping CLIP-seq data. You need to set up the alignment tool in the `config/user_config.yaml` file correctly. It is **recommended** to specify the number of threads in `config/user_config.yaml` file by adding `threads_mapping: N`, or you can simply add `-j N` parameter in the CLIPipe command. The other detial parameters for alignment can be found in `config/default_config.yaml`.

```bash
cd /home/CLIPipe_user/clipipe/clipipe_demo/general/;
clipipe -u ./config/user_config.yaml mapping
```

```text
Note:
    The output folder `output_human_hg38/mapping_bowtie/` contains alignment results using bowtie.
    The output folder `output_human_hg38/mapping_bwa/` contains alignment results using bwa.
    The output folders `output_human_hg38/mapping_novoalign/` contain alignment results using novoalign.
```

### Peak calling

CLIPipe provides multiple peak calling methods for identifying recurring fragments of CLIP-seq data.

<table>
    <tr>
        <th></th>
        <th>Method-specific</th>
        <th>Non-specific</th>
    </tr>
    <tr>
        <td>HITS-CLIP</td>
        <td>CTK</td>
        <td rowspan="15">Piranha</td>
    </tr>
    <tr>
        <td>PAR-CLIP</td>
        <td>PARA suite</td>
    </tr>
    <tr>
        <td>iCLIP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>iCLAP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>eCLIP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>4sU-iCLIP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>urea-iCLIP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>BrdU-CLIP</td>
        <td>CTK, PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>Fr-iCLIP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>FAST-iCLIP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>irCLIP</td>
        <td>PureCLIP, iCLIPro,  iCount</td>
    </tr>
    <tr>
        <td>seCLIP</td>
        <td>PureCLIP</td>
    </tr>
    <tr>
        <td>uvCLAP</td>
        <td>JAMM</td>
    </tr>
    <tr>
        <td>FLASH</td>
        <td>PureCLIP</td>
    </tr>
    <tr>
        <td>dCLIP</td>
        <td>PeakRanger</td>
    </tr>
</table>


```bash
cd /home/CLIPipe_user/clipipe/clipipe_demo/general/;
clipipe -u ./config/user_config.yaml peak_calling;    # Please choose from Piranha(mapping method: biwtie) CTK(mapping method: novoalign) PureCLIP(mapping method: biwtie) parclip_suite(do not need mapping step)
```

```text
Note:
    The output folders `output_human_hg38/peak_calling_piranha/` contain alignment results using piranha.
    The output folder `output_human_hg38/peak_calling_CTK/` contains peak calling results using CTK.
    The output folders `output_human_hg38/peak_calling_pureclip/` contain alignment results using pureclip.
    The output folders `output_human_hg38/peak_calling_parclip_suite/` contain alignment results using parclip_suite.
```

Several other peak calling tools can be used in the CLIPipe docker directily:

```bash
# iCLIPro
$ iCLIPro [options] in.bam

# iCount
$ iCount [-h] [-v] ...

# JAMM
$ JAMM.sh --help

# PeakRanger
$ peakranger <command> <arguments>

# clipcontext
$ clipcontext [-h] [-v] {g2t,t2g,lst,int,exb,eir} ...
```

### Motif discovery

The motif discovery function can be used directly in the CLIPipe docker:

For HOMER, the demo script like this:

```bash
# input: ${sample_id}.all_peak.bed

# 1. split training and test dataset
perl /home/CLIPipe_user/clipipe2/clipipe_software/bin/homer/1.split.pl ${sample_id}.all_peak.bed;

# 2. prepare training and test fasta
perl /home/CLIPipe_user/clipipe2/clipipe_software/bin/homer/2.prepare_Homer.pl ${sample_id} training ${genome_fasta};
perl /home/CLIPipe_user/clipipe2/clipipe_software/bin/homer/2.prepare_Homer.pl ${sample_id} test ${genome_fasta};

# 3. Run Homer on training dataset
findMotifs.pl ${sample_id}.training_peak.fa fasta Homer_training_output -len 4,5,6,7,8,9,10 -rna    # the number of len could change

# 4. Run Homer on test dataset
mkdir Homer_test_output;
findMotifs.pl ${sample_id}.test_peak.fa fasta Homer_test_output -rna -find Homer_training_output/homerMotifs.all.motifs > Homer_test_output/count.txt;

```

For MEME, the demo script like this:

```bash
# input: ${sample_id}.all_peak.bed

# 1. split training and test dataset
perl /home/CLIPipe_user/clipipe2/clipipe_software/bin/meme/1.split.pl ${sample_id}.all_peak.bed;

# 2. prepare training and test fasta
perl /home/CLIPipe_user/clipipe2/clipipe_software/bin/meme/2.prepare_MEME.pl ${sample_id} training ${genome_fasta};
perl /home/CLIPipe_user/clipipe2/clipipe_software/bin/meme/2.prepare_MEME.pl ${sample_id} test ${genome_fasta};

# 3. run MEME on training dataset
mkdir MEME_output
meme ${sample_id}.training_peak.fa -o MEME_output -dna -minw 4 -maxw 10 -nmotifs 25; # the number of minw, maxw and nmotifs could change

# 4. run FIMO on test dataset
cat MEME_output/meme.txt | sed 's/10.0e+000/1.0e+001/g' | sed 's/10.0e+001/1.0e+002/g' | sed 's/10.0e+002/1.0e+003/g' | sed 's/10.0e+003/1.0e+004/g' | sed 's/10.0e+004/1.0e+005/g' | sed 's/10.0e+005/1.0e+006/g' | sed 's/10.0e+006/1.0e+007/g' | sed 's/10.0e+007/1.0e+008/g' | sed 's/10.0e+008/1.0e+009/g' | sed 's/10.0e+009/1.0e+010/g' | sed 's/10.0e+010/1.0e+011/g' > meme.txt;
fimo --thresh 0.01 -o FIMO_output meme.txt ${sample_id}.test_peak.fa; # the number of thresh could change

```

Other related tools are also provided:

```bash
# PhyloGibbs
$ phylogibbs-mp [-m motifwidth] input_seqfile [input_seqfile2 ...]

# STREME
$ streme [options]
```

## Copyright and License Information

Copyright (C) 2021 Tsinghua University, Beijing, China

This program is licensed with commercial restriction use license. Please see the [LICENSE](https://github.com/ShangZhang/clipipe/blob/main/LICENSE) file for details.
