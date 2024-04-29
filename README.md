# Somatic copy number alteration (SCNA) pattern identification in shallow whole-genome sequencing data

## Introduction

**dincalcilab/scnapattern** is a bioinformatics pipeline that assigns somatic copy number alterations (SCNA) patterns (S, stable; U, unstable; HU, highly unstable) as described in [Copy number alterations in stage I epithelial ovarian cancer highlight three genomic patterns associated with prognosis](https://doi.org/10.1016/j.ejca.2022.05.005) by Pesenti, Beltrame, *et al*.

This pipeline requires shallow whole-genome sequencing (sWGS) data with absolute copy number and ploidy estimates. Currently segmentation outputs from [ASCAT.sc](https://github.com/VanLoo-lab/ASCAT.sc), [ACE](https://www.bioconductor.org/packages/release/bioc/html/ACE.html), and [ichorCNA](https://github.com/GavinHaLab/ichorCNA/) are supported natively.

The pipeline includes these steps:

1. Calculate SCNA patterns for given inputs
2. Output a table with patterns and parameters for each analyzed sample
3. Generate a QC summary ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,filename,ploidy,format
Sample_1,Sample1.cna.seg,2,ichorcna
```

Each row represents a sample name, the associated absolute path to the segment file, the ploidy of the sample, and the data format (either `ascat`, `ace`, or `ichorcna`).


Now, you can run the pipeline using:

```bash
nextflow run dincalcilab/scnapattern \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --genome <GENOME> \
   --genomestyle <GENOMESTYLE>
```

Where `GENOME` is either `hg19` or `hg38` and `genomestyle` is either `ucsc` (`chr` prefix for chomosomes) or `ncbi` (no `chr` prefix). Note that you **must** use a Nextflow configuration profile that supports Docker or Singularity images, as some tools are only provided by containers.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

dincalcilab/scnapattern was originally written by Luca Beltrame (@lbeltrame).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use dincalcilab/scnapattern for your analysis, please cite it using the following doi: [10.1016/j.ejca.2022.05.005](https://doi.org/10.1016/j.ejca.2022.05.005)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
