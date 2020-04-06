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

### Public Docker images
If you would like an additional Docker image, please first check [if it is in the Biocontainers registry](https://biocontainers.pro/#/registry).
If it's not, it may be in the `hpobiolab` dockerhub. All tools in the `dockerfiles` directory are available via `docker pull hpobiolab/<toolname>`, for
example.

The public docker images used for tools should come from the hpobiolab (or another vetted, public repo) to ensure both security and reliability.

Docker images currently being used are here: [https://hub.docker.com/orgs/hpobiolab/repositories](https://hub.docker.com/orgs/hpobiolab/repositories)


### Running workflows
WDL requires cromwell to run locally.

To download cromwell: wget https://github.com/broadinstitute/cromwell/releases/download/49/cromwell-49.jar

To download womtool: wget https://github.com/broadinstitute/cromwell/releases/download/49/womtool-49.jar

To run in the cloud, you can run in Broad's Terra environment or using Google Cloud Project directly. Instructions for doing so will follow shortly.

To run Nextflow pipelines, please see the `nextflow` directory.


## Original proposal
The pangenomics channel is working on generating assembly-based pangenomes of SARSCov2 genomes. Since we already have a reference genome (including a GFF file of ORF annotations), I thought it might be useful to build analysis pipeline(s) that can operate in parallel or downstream of the assembly pangenome.

NextStrain already does things like convert the RNA/cDNA sequences to amino acids. I was thinking we could use either their tooling or our own to produce some automatically-generated reports of variable sites on the genome / proteome. We can also provide these annotations as GFA paths to incorporate into the pangenome, facilitate read alignment to ref genome / pangenome, or filter reads against viral or host references using Kraken / rkmh.

I'm most comfortable in WDL (which runs in Broad's Terra, DNANexus via dxWDL, and using Google's Pipelines API), but we could use any of the workflow languages in reality. I think this would be a good project for folks wanting to work in shell, WDl, python, docker, and certainly R as well.

Scope-wise, it's probably best to start with a single workflow that annotates variable sites, then try to build one that aligns reads and reports whether a new strain has novel variation at these (or other) sites. Filtering workflows could be a component of this workflow.

## Pipeline(s)

### Read filtering (in progress)

### Read-to-reference alignment (in progress)

### Variable site detection

### Pangenome generation with minimap2, seqwish, and odgi (in progress)
We've wrapped the commands to build the odgi-based SARS-CoV2 pangenome and
visualization, described [here](https://github.com/virtual-biohackathons/covid-19-bh20/wiki/Pangenome), into a WDL description that can be run via Cromwell.

If you have docker installed locally, the workflow should run after pulling the proper images
(all built from dockerfiles in the `dockerfiles` directory).

To run locally, you'll need to make sure Cromwell is downloaded. There's a directory for this but we don't include it in the repo as the JAR file is rather large.

Next, you'll need to create an `inputs.json` file that tells cromwell where your reads are.
For a FASTA file of sequences named "seqs.fa" in the current directory,
your file would look like so:
```
{
   PangenomeGenerate.inputReads="seqs.fa"
}
```

Save this file as `inputs.json`. Then, you can run cromwell locally like so:  

```
java -jar cromwell/cromwell-49.jar run -i inputs.json workflow/pangenome-generate.wdl
```


If you're running in the cloud,
you'll want to set up a google cloud project / billing account
and fill in the missing fields in the
EXAMPLE conf file. Steps roughly outlined below:
1. Open EXAMPLE.CROMWELL.PAPI.conf in your favorite text editor
2. Replace the `<your-google-project-here>` placeholders (there should be three):
```
.
.
.
engine {
  filesystems {
    gcs {
      auth = "application-default"
      project = "<YOUR PROJECT ID HERE>"
    }
  }
}
.
.
.
config {
        // Google project
        project = "<YOUR PROJECT ID HERE>"

        // Base bucket for workflow executions
.
.
.
filesystems {
          gcs {
            // A reference to a potentially different auth for manipulating files via engine functions.
            auth = "application-default"
            project = "<YOUR PROJECT ID HERE>"
          }
        }
.
.
.
```
3. Replace the Google bucket ID place holder with your google bucket name:
```
.
.
.
      project = "<YOUR PROJECT ID HERE>"

        // Base bucket for workflow executions
        root = "gs://<YOUR GOOGLE STORAGE BUCKET NAME HERE>/cromwell-execution"
.
.
.
```


Next, you'll also need to copy your input reads to your cromwell google cloud bucket:
```
gsutil cp seqs.fa gs://<your-google-bucket>/
```

and rather than use the file name in the inputs.json,
put the absolute path in your cromwell bucket:
```
{
    PangenomeGenerate.inputReads="gs://<your-google-bucket>/seqs.fa"
}
```

Once you've done that, you can run like so:

```
java -Dconfig.file=cromwell/EXAMPLE.CROMWELL.PAPI.conf -jar cromwell.jar -i inputs.json workflow/pangenome-generate.wdl
```

Cromwell will run in the foreground and report status to the terminal. When it's done, you should recieve a "Success" message, and your files will be in your google bucket.


### Read-to-graph alignment with vg (in progress)


## Examples

## Building custom Docker images

The Makefile provides some shorthand for building/pushing dockerfiles. From the main project directoy,
a `make build APP=<APPNAME>` will build an app named "APPNAME" from a dockerfile named "APPNAME.Dockerfile."

```
make build APP=<dockerfile basename>

## Example for samtools.Dockerfile
make build APP=samtools

## To push to a public repo:
make push APP=samtools
```

Here's a concrete example, based on the BWA dockerfile in the dockerfiles directory.
```
make build APP=bwa
make push APP=bwa
```

The resulting docker image is then pushed to [https://hub.docker.com/repository/docker/hpobiolab/bwa](https://hub.docker.com/repository/docker/hpobiolab/bwa)

## Reference genomes and annotations
The `refs` directory contains the SARS-COV2 reference genome (from GenBank / RefSeq) as well
as the GFF file of genomic features. In addition, a pre-built BWA index set is included for read mapping.

## Data
The data directory contains links for downloading raw read data from SRA. Raw reads should not be included in the repo as they tend to break git.


