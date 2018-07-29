# auto-asm

Automated assembly pipeline for PacBio long read datasets.


## Setup

### Install miniconda

Be sure to accept the option to automatically modify your `.bashrc` or manually
update it so that your installed miniconda packages are in your `PATH`.

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

### Install smrtcmds

Download and install SMRT Tools (smrtcmds) from
[PacBio](https://www.pacb.com/support/software-downloads/). SMRT Link is not
required; the installer offers a SMRT Tools-only installation option. **Do not**
forget to put the path to the `smrtcmds/bin` folder in your configuration.

auto-asm has been tested with PacBio SMRT Tools 5.1.

### Clone and configure auto-asm

Clone the repository and copy the sample configuration file to a new file.
A sensible location would be the root of the desired auto-asm working/output
directory.

**Note:** auto-asm was originally developed to assemble _A. thaliana_ genomes
and the default resource configuration (cores/memory) is set with this in mind.
If you are assembling significantly smaller or larger genomes, you will have
to adjust these values. See the [Configuration Guide](#configuration-guide)
below.

```
$ git clone https://github.com/weigelworld/auto-asm
$ cd auto-asm
$ cp config.yaml.sample config.yaml
$ vim config.yaml
```

### Install snakemake

```
$ conda env create --name "auto-asm" --file envs/auto-asm.yaml
```

### Run auto-asm

auto-asm comes with a sample bash script (`auto-asm.sh.sample`) for
running the pipeline on an SGE cluster. You may need to edit the cluster
commands in the sample script to work in your cluster environment. Also review
`cluster.json` and make sure that it is appropriate for your cluster setup,
especially the cluster output and error file naming patterns.

Copy the `auto-asm.sh.sample` script and edit the working directory and
`config.yaml` location as needed. Run the bash script to start the auto-asm
pipeline.

```
$ cp auto-asm.sh.sample auto-asm.sh
$ vim auto-asm.sh
$ chmod +x auto-asm.sh
$ source activate auto-asm
$ ./auto-asm.sh
```


## Configuration Guide

If you copied the sample configuration file, the required rule parameters are
already set to their defaults. You should only need to replace the assemblies
and smrtcmds\_bin path but you can also change the resources and other
parameters to suit your usecase. The sample configuration is a good guideline
for how to format the configuration in YAML, for those unfamiliar with it.

For users that understand JSON schemas, the configuration schema is located at
`schemas/config.schema.yaml`. Configuration files are validated against this
schema at runtime.

The 'assemblies' property is a list of mappings between alphanumeric assembly
ids and their configuration objects. Each assembly configuration must, at
minimum, have `genome_size`, `eukaryotic`, and `long_read_paths` defined.

Assembly configuration options:

- name
  - optional
  - nicely-formatted assembly name (e.g. Col-0)
- reference
  - optional
  - path to reference genome in FASTA format
- genome\_size
  - **required**
  - genome size in bp, with optional k, m, or g suffixes (e.g. 120m)
- eukaryotic
  - **required**
  - whether the organism is eukaryotic or not
- long\_read\_paths
  - **required**
  - list of paths to PacBio long read BAM files
- paired\_short\_read\_paths
  - optional
  - list of pairs of paths to FASTQ files of a paired-end short read dataset

