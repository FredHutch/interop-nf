# InterOp QC Reporting - Nextflow

Workflow used to generate QC metrics from the output of an Illumina sequencing run

## Invoking the Workflow

The workflow can be run on the contents of an Illumina run folder with the following
command:

```#!/bin/bash

runDirectory=<path to run directory>
reportsDirectory=<path to report output folder>

nextflow run FredHutch/InterOp-nf --runDir ${runDirectory} --reportsDir ${reportsDirectory}

```

The run folder must contain the `InterOp/` folder and the `RunInfo.xml` and `RunParameters.xml` files.

## Other information

Docker image used:
https://github.com/FredHutch/docker-python-interop