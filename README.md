# InterOp QC Reporting - Nextflow

Workflow used to generate QC metrics from the output of an Illumina sequencing run.

[![Run tests](https://github.com/FredHutch/interop-nf/actions/workflows/tests.yml/badge.svg)](https://github.com/FredHutch/interop-nf/actions/workflows/tests.yml)

[![codecov](https://codecov.io/gh/FredHutch/interop-nf/branch/main/graph/badge.svg?token=I5RTB5TXLA)](https://codecov.io/gh/FredHutch/interop-nf)

## Invoking the Workflow

The workflow can be run on the contents of an Illumina run folder with the following
command:

```#!/bin/bash

runDirectory=<path to run directory>
reportsDirectory=<path to report output folder>

nextflow run FredHutch/InterOp-nf --runDir ${runDirectory} --reportsDir ${reportsDirectory}

```

The run folder must contain the `InterOp/` folder and the `RunInfo.xml` and `RunParameters.xml` files.

## Workflow output

[View more info here](sample_output/README.md).

## Other information

Docker image used:
https://github.com/FredHutch/docker-python-interop