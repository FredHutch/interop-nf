#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Initialize command line flag parameters
params.runDir = false
params.reportsDir = false
params.help = false

// Command to print the help message
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/InterOp-nf --runDir <> --reportsDir <>
    
    Required Arguments:
      --runDir               Path to illumina run directory
      --reportsDir           Path to the report output

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
// or does not provide both --runDir and --reportsDir
if (params.help || params.runDir == false || params.reportsDir == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// DEFINE PROCESSES
process extractStats {
    container "quay.io/hdc-workflows/python-interop"
    label 'io_limited'
    errorStrategy 'finish'
    publishDir "${params.reportsDir}", mode: 'copy', overwrite: true
    
    input:
    path xml_inputs, stageAs: 'runDir/'
    path interop_inputs, stageAs: 'runDir/InterOp/'

    output:
    file "run_metrics.json"
    file "run_metrics.html"
    file "percent_base.pdf"
    file "max_intensity.pdf"
    file "occupancy.pdf" optional true

    """#!/bin/bash

set -Eeuo pipefail

run_summary.py runDir/

plot_percent_base.py runDir/

plot_tile_intensity.py runDir/

plot_occupancy.py runDir/
    """

}

// DEFINE THE WORKFLOW
workflow {

    // Extract the JSON with quality metrics from the
    // xml files in the input folder
    extractStats(
        Channel.fromPath("${params.runDir}/*xml").toSortedList(),
        Channel.fromPath("${params.runDir}/InterOp/*bin").toSortedList(),
    )

}