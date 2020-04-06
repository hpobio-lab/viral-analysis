nextflow
---

To install Nextflow, please see https://www.nextflow.io/index.html#GetStarted.


To run locally, all prerequistes must be installed

```bash
$ nextflow run pangenome-generate.nf
```

otherwise run via Docker

```bash
$ nextflow run pangenome-generate.nf -profile docker
```

By default these workflows look for data in `../data` and copy results to an
output directory based on the simple name of input file (e.g. `SARS-CoV-2/` for
`SARS-CoV-2.genbank.20200329.complete.fasta`).

To run on other executors (such as SGE, SLURM, Kubernetes, AWS Batch, and more),
please see https://www.nextflow.io/docs/latest/executor.html.
