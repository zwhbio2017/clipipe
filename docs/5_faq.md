## Frequently Asked Questions

### How to install and set up the environment?

We have packaged everything into CLIPipe docker image, you can [download the raw image](http://clipipe.ncrnalab.org/CLIPipe_v1.0.3_.tar.gz) or get it from [docker hub](https://hub.docker.com/repository/docker/shangzhang/clipipe)

### How to use CLIPipe?

CLIPipe is a comprehensive quality control and analysis pipeline for CLIP-seq data. We use snakemake for parallel running and further integrate snakemake pipeline into one single command.

The basic usage of CLIPipe are described [here](3_basic_usage.md). Basically you should complete the following steps before running the command:

-   Set up environment and install docker
-   Prepare genome and annotation
-   Install CLIPipe
-   prepare input files in right file path
-   set up configuration

Then you can run the command, you can specify the module you want to run and dataset you provide.

```bash
clipipe ${step_name} -d ${dataset}

# Note:
#   ${step_name} is one of the step listed in 'positional arguments'.
#   ${dataset} is the name of your dataset that should match the prefix of your configuration file described in the following section.
```

### What is Snakemake?

The [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system is a tool to create reproducible and scalable data analyses. We have hide the details of snakemake and you only need to run one single command. However you can customize some of the codes if you are familiar with snakemake.

### How to set configurations in config file

There are many parameters to be specified. You'd better making a new copy of the config file in config directory. For example you can make one copy of `user_config.yaml`. Then rename the file to `config/${dataset}.yaml`.

Other parameters are defined in `default\_config.yaml`. You may also change parameters.

### Why PureCLIP run so slowly and often use up the memory?

It is possibly because PureCLIP builds models for all chromosomes in the genome fasta files, regardless of whether the bam file contains alignment on the chromosome. Thus, PureCLIP would build models for many chromosomes and scaffolds that even have no reads on it, which would be a waste of time and memory. If you encounter such problems, consider specifying the "-iv" parameter when running PureCLIP in such format: `-iv 'chr1;chr2;chr3'` (must match the chromosome name in the genome fasta file). In this way, PureCLIP would only build models on the 3 defined chromosomes, thus making it faster to process the data.

### When bugs appear?

The quickest way is to create a [new issue](https://github.com/ShangZhang/clipipe/issues)

If you want us to add more functions in CLIPipe, please create a [new issue](https://github.com/ShangZhang/clipipe/issues)

### Why did some jobs occasionally fail with no extra error message except 'CalledProcessError'?

The most possible cause is no available memory. You can confirm the problem by running the linux command `dmesg` or open and examine the most recent message. Message like "Out of memory: Kill process \*\* or sacrifice child" clearly indicates that memory problem occurred. Some jobs \(e.g. mapping using bwa, samtools sort, bedtools sort\)requires large amount memory especially when the input number of reads is large. Try to reduce the number of parallels with the `-j` option or set memory limit for a particular job if you run the jobs on a computer cluster.

### How to rerun downstream steps from a specific step in a pipeline?

Sometimes we need to rerun a pipeline from a step, usually after changing the configuration file. Snakemake is not aware of changes in configuration file and we need to rerun the pipeline by ourselves. The `--forcerun` option in snakemake allows rerunning a step and all steps that depend on the output files of the step. For example, to rerun the `peak_calling` step, just run:

```bash
clipipe peak_calling -d $dataset --forcerun
```
