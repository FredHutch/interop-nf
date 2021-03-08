#!/usr/bin/env python3

import logging
import os
import sys

import matplotlib.pyplot as plt
from interop import py_interop_run_metrics, py_interop_run, py_interop_plot

logger = logging.getLogger(__name__)


def plot_percent_base(run_folder: str, output_pdf="percent_base.pdf"):
    """
    Plots the base % across each cycle. Each line represents a different base.
    Reference lines are added for each read.

    Base %: The percentage of clusters for which the selected base (A, C, T, or G) has been called.
    """
    # Initialize interop objects
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)

    # Read from the run folder
    run_metrics.read(run_folder, valid_to_load)

    logger.info('Generating % base plot')
    plot_data = py_interop_plot.candle_stick_plot_data()
    options = py_interop_plot.filter_options(run_metrics.run_info().flowcell().naming_method())
    py_interop_plot.plot_by_cycle(run_metrics, "BasePercent", options, plot_data)

    # Plot each base
    for base_index in range(plot_data.size()):
        line_data = plot_data.at(base_index)
        x = [line_data.at(i).x() for i in range(line_data.size())]
        y = [line_data.at(i).y() for i in range(line_data.size())]
        plt.plot(x, y, color=line_data.color(), linewidth=0.5, label=line_data.title())

    # Plot reference lines for reads
    read_vector = run_metrics.run_info().reads()
    for read_index in range(read_vector.size()):
        read_name = f'R{read_vector[read_index].number()}'
        cycle_start = read_vector[read_index].first_cycle()
        plt.axvline(x=cycle_start, color='purple', linestyle='--', linewidth=0.35)
        plt.text(cycle_start, plt.gca().get_ylim()[1], read_name, fontsize=8, color='purple')

    # Plot settings
    axes_data = plot_data.xyaxes()
    plt.xlabel(axes_data.x().label(), fontsize=10)
    plt.ylabel(axes_data.y().label(), fontsize=10)
    plt.title(plot_data.title(), fontsize=10)
    plt.legend()
    plt.ylim([axes_data.y().min(), axes_data.y().max()])
    plt.xlim([axes_data.x().min(), axes_data.x().max()])

    # Save figure
    plt.savefig(output_pdf)


# If this file is being run as a standalone script
if __name__ == "__main__":

    # Get the input file path from the first argument
    assert len(sys.argv) > 1, "Please provide input path"
    input_path = sys.argv[1]
    assert os.path.exists(input_path), f"Input path must exist ({input_path})"

    plot_percent_base(input_path)
