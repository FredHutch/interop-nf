# InterOp QC Reporting - Nextflow

Workflow used to generate QC metrics from the output of an Illumina sequencing run

## Invoking the Workflow

The workflow can be run on the contents of an `InterOp` folder with the following
command:

```#!/bin/bash

INPUT=<path to input InterOp folder>
OUTPUT=<path to output folder>

nextflow run FredHutch/InterOp-nf --input ${INPUT} --output ${OUTPUT}

```
