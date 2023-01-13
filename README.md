# AlignerCompare

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/alignercompare)

## Introduction

Compare
## Pipeline summary

1. Create indicies, process annotation files, create splice junctions for
the various aligners if not provided
2. Align with various aligners. Currently set up:
  - minimap2
  - STAR
  - hisat2
3. QC the aligners
  - samtools stats,idxstats,flagstat
  - RSEQC
  - Subreads/featurecounts
4. Collate QC reports with MultiQC

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Pull this repo

4. Run the test. I suggest creating a file with the following cmd:

  ```bash
  #!/bin/bash

  export NXF_SINGULARITY_CACHEDIR=singularity

  nextflow run aligner_compare/main.nf \
      -profile test,singularity
  ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

To do this, you will need to create a sample sheet `csv`. Here is an example of
a samplesheet to run HG002 through:

```raw
sample,fastq_1,fastq_2
hg002,/home/oguzkhan/projects/testing_aligncompare/NA24385.filtered.fq.gz,
```

Adjust the path according to your system. Note that the SQANTI fasta is converted
to a fastq file -- I have a script set up in the oop_ify branch of isocomp to
do this.

Next, set up the input parameters json file:

```json
{
   "input": "samplesheet.csv",
   "outdir": "results",
   "fasta": "/home/oguzkhan/ref/human/GRCh38.primary_assembly.genome.fa.gz",
   "gtf": "/mnt/isocomp/genome/gencode/release_42/gencode.v42.annotation.gtf"
}

```

And the aligner configuration file. I call this `aligner_settings.config` but
the name doesn't matter:

```raw
// minimap2
// note that -a is hard coded into the bam output cmd of the
// nextflow module minimap command. no need to include -a here
process {

    withName: MINIMAP2_ALIGN {
        ext.args = [
            '-x splice:hq',
            '-uf'
        ].join(' ').trim()
    }

}

// hisat2
process {

   withName: HISAT2_ALIGN {
       ext.args = [
           '--min-intronlen 20',
           '--max-intronlen 1000000'
       ].join(' ').trim()

   }

}

// STAR
process {

    withName: STAR_ALIGN {
        ext.args = [
            '--readFilesCommand zcat',
            '--STARlong',
            '--twopassMode Basic',
            '--outFilterType BySJout',
            '--outFilterMultimapNmax 20',
            '--alignSJoverhangMin 8',
            '--alignSJDBoverhangMin 1',
            '--outFilterMismatchNmax 999',
            '--outFilterMismatchNoverReadLmax 0.04',
            '--alignIntronMin 20',
            '--alignIntronMax 1000000',
            '--alignMatesGapMax 1000000'
        ].join('').trim()
    }

}

```

A run command for this would look like so -- I paste this into a file:

```bash
#!/bin/bash

export NXF_SINGULARITY_CACHEDIR=singularity

# path to the aligner_compare/main.nf file -- make sure correct for your system
# you may want to change the profile to docker or whatever other dependency
# manager you're using
nextflow run aligner_compare/main.nf \
    -profile singularity \
    -c aligner_settings.config \
    -params-file params.json \
    --max_memory 20.GB \
    --max_cpus 10 \
    -resume

```

## Credits

nf-core/alignercompare was originally written by chase mateusiak.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/alignercompare for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
