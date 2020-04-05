Viral analysis workflows in the cloud (part of COVID-19 Biohackathon)
--------------------
April 2020

## Overview
This repo holds workflows for analyzing viral genomes. It is
particularly tailored to SARS-CoV2 as part of the [COVID-19 Biohackathon](https://github.com/virtual-biohackathons/covid-19-bh20).

To chat with devs, visit the biohackathon slack and see the Workflows channel.

## Installation / requirements

### Development
You'll need Docker installed locally do build dockerfiles, as well as GNU make to use the makefile.


### Running workflows
WDL requires cromwell to run locally.

To download cromwell: wget https://github.com/broadinstitute/cromwell/releases/download/49/cromwell-49.jar

To download womtool: wget https://github.com/broadinstitute/cromwell/releases/download/49/womtool-49.jar

To run in the cloud, you can run in Broad's Terra environment or using Google Cloud Project directly. Instructions for doing so will follow shortly.



## Original proposal
The pangenomics channel is working on generating assembly-based pangenomes of SARSCov2 genomes. Since we already have a reference genome (including a GFF file of ORF annotations), I thought it might be useful to build analysis pipeline(s) that can operate in parallel or downstream of the assembly pangenome.

NextStrain already does things like convert the RNA/cDNA sequences to amino acids. I was thinking we could use either their tooling or our own to produce some automatically-generated reports of variable sites on the genome / proteome. We can also provide these annotations as GFA paths to incorporate into the pangenome, facilitate read alignment to ref genome / pangenome, or filter reads against viral or host references using Kraken / rkmh.

I'm most comfortable in WDL (which runs in Broad's Terra, DNANexus via dxWDL, and using Google's Pipelines API), but we could use any of the workflow languages in reality. I think this would be a good project for folks wanting to work in shell, WDl, python, docker, and certainly R as well.

Scope-wise, it's probably best to start with a single workflow that annotates variable sites, then try to build one that aligns reads and reports whether a new strain has novel variation at these (or other) sites. Filtering workflows could be a component of this workflow.

## Pipeline(s)

## Examples

## Building custom Docker images
To build a docker file in the dockerfiles directory, try:
```
make build APP=<dockerfile basename>

## Example for samtools.Dockerfile
make build APP=samtools

## To push to a public repo:
make push APP=samtools
```

## Reference genomes and annotations

## Data


