#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Initialize command line flag parameters
params.input = false
params.output = false
params.help = false

// Command to print the help message
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/InterOp-nf --input <> --output <>
    
    Required Arguments:
      --input               Path to input folder
      --output              Path to output folder

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
// or does not provide both --input and --output
if (params.help || params.input == false || params.output == false){
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
    publishDir "${params.output}", mode: 'copy', overwrite: true
    
    input:
    path xml_inputs, stageAs: 'input/'
    path interop_inputs, stageAs: 'input/InterOp/'

    output:
    file "qc_metrics.json"
    file "percent_base_plot.png"
    file "run_stats.html"
    file "max_intensity.pdf"

    """#!/bin/bash

set -Eeuo pipefail

interop_helper.py input/

plot_tile_intensity.py input/
    """

}

// DEFINE THE WORKFLOW
workflow {

    // Extract the JSON with quality metrics from the
    // xml files in the input folder
    extractStats(
        Channel.fromPath("${params.input}**xml").toSortedList(),
        Channel.fromPath("${params.input}**InterOp/*bin").toSortedList(),
    )

}