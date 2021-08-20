[TOC]
# CLIPipe

CLIPipe(CLIP-seq Pipeline) is an integrated pipeline for analyzing CLIP sequencing data.

![Pipeline of Tutorial](tutorial/CLIPipe_pipeline.png)

The CLIPipe workflow consists of:

- Pre-processing function:
  - Quality control of the raw CLIP-seq reads.
  - Remove adapter of the sequencing data.
  - Filter low quality reads of the sequencing data.
  - Collpase duplicates of the sequencing data.
  - Remove barcode of the sequencing data.

- Alignment function :
  - mapping sequencing data to reference genome using bowtie, bwa and novoalign

- Peak calling function :
  - Peak calling of the CLIP-seq data using Piranha, CTK, PureCLIP, iCLIPro, iCount, JAMM, PEAKachu, PeakRanger and clipcontext

- Motif discovery function :
  - Motif discovery of the CLIP-seq data using HOMER, PhyloGibbs, MEME, GraphProt, DREME and STREME

## Tutorial

<iframe height='265' scrolling='no' title='CLIPipe Tutorial' src='https://clipipe.readthedocs.io/en/latest/' frameborder='no' allowtransparency='true' allowfullscreen='true' style='width: 100%;'></iframe>


## Requirements
### software
* docker

## Installation

  1. All required software and packages are already installed in docker, so there are no more requirements. You can execute to get the docker `CLIPipe_v3` container:
  
    ```bash
    docker run --name=CLIPipe_v3 -t -d -h CLIPipe_docker --restart unless-stopped -v /lulab/lustre2/zhangshang/work/software/clipipe:/home/CLIPipe_user/clipipe zs/clipipe:v1 /bin/bash
    ```

2. To show the docker `CLIPipe_v3` container, you can execute:

   ```bash
   docker container ls
   ```

3. To execute the `CLIPipe_v3` container, you can execute:

   ```bash
   docker exec -it CLIPipe_v3 bash
   ```

4. After entering the container, please change the user to `CLIPipe_user`

   ```bash
   su CLIPipe_user;
   cd ~
   ```

5. To test the installation and get information about the command-line interface of CLIPipe, you can execute:

   ```bash
   clipipe --help
   ```

   A helper message is shown like this:

   ```
   usage: clipipe [-h] --user_config_file USER_CONFIG_FILE
                  {pre_process,mapping,peak_calling}
   
   CLIPipe: A comprehensive quality control and analysis pipeline for CLIP highthroughput sequencing data
   ==================================================================================
   CLIPipe is a Python module and associated command-line interface (CLI), which provides all the commands needed to process protein-RNA CLIP interaction data and to identify and quantify sites of protein-RNA interactions on RNA.
   
   CLIPipe's main input are FASTQ files with iCLIP sequencing data, its main output are BED files with identified and quantified cross-linked sites.
   
   A number of analyses are included in CLIPipe that provide insights into the properties of
   protein-RNA interaction.
   
   positional arguments:
     {pre_process,mapping,peak_calling}
   
   optional arguments:
     -h, --help            show this help message and exit
     --user_config_file USER_CONFIG_FILE, -u USER_CONFIG_FILE
                           the user config file
   
   ==================================================================================
   For additional help or support, please visit https://github.com/ShangZhang/clipipe
   ```

## Basic Usage
The basic usage of CLIPipe is:
```bash
clipipe ${step_name} -d ${dataset}
```

> **Note:**
>
> - `${step_name}` is one of the step listed in 'positional arguments'.
> - `${dataset}` is the name of your dataset that should match the prefix of your configuration file described in the following section.


## Usage

### Pre-processing



## Copyright and License Information

Copyright (C) 2021 Tsinghua University, Beijing, China

This program is licensed with commercial restriction use license. Please see the [LICENSE](https://github.com/ShangZhang/clipipe/blob/main/LICENSE) file for details.

