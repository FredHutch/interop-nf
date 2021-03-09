#!/usr/bin/env python3

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from interop import py_interop_run_metrics, py_interop_run, py_interop_table
from matplotlib.backends.backend_pdf import PdfPages


def plot_occupancy(run_folder, output_pdf="occupancy.pdf"):

    # Initialize interop objects
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid_to_load[py_interop_run.ExtendedTile] = 1
    valid_to_load[py_interop_run.Tile] = 1
    valid_to_load[py_interop_run.Extraction] = 1

    # Read from the run folder
    run_metrics.read(run_folder, valid_to_load)

    # Create the columns
    columns = py_interop_table.imaging_column_vector()
    py_interop_table.create_imaging_table_columns(run_metrics, columns)

    headers = []
    for i in range(columns.size()):
        column = columns[i]
        if column.has_children():
            headers.extend(
                [column.name()+"("+subname+")" for subname in column.subcolumns()])
        else:
            headers.append(column.name())

    column_count = py_interop_table.count_table_columns(columns)
    row_offsets = py_interop_table.map_id_offset()
    py_interop_table.count_table_rows(run_metrics, row_offsets)
    data = np.zeros((row_offsets.size(), column_count), dtype=np.float32)
    py_interop_table.populate_imaging_table_data(
        run_metrics, columns, row_offsets, data.ravel()
    )

    # Make a DataFrame
    df = pd.DataFrame(data, columns=headers)

    # If there is no data
    if df.shape[0] == 0:
        # Stop
        return

    # % Occupied only available on NovaSeq
    if "% Occupied" not in df:
        print("Occupancy plot skipped, no data available")
        return

    # Otherwise

    # Make a PDF
    with PdfPages(output_pdf) as pdf:

        # Make a few different types of plots
        for x, y, hue in [
            ("% Occupied", "% Pass Filter", "Tile"),
            ("% Occupied", "% Pass Filter", "Lane"),
            ("% Occupied", "% Pass Filter", "Cycle"),
            ("% Occupied", "% Pass Filter", "Read"),
        ]:
            g = sns.scatterplot(
                data=df,
                x=x,
                y=y,
                hue=hue,
                alpha=0.5,
                linewidth=0,
            )
            plt.xlim([0, 100])
            plt.ylim([0, 100])
            plt.legend(title=hue, bbox_to_anchor=[1.2, 0.9])
            plt.tight_layout()
            pdf.savefig(bbox_inches="tight")
            plt.close()

# If this file is being run as a standalone script
if __name__ == "__main__":

    # Get the input file path from the first argument
    assert len(sys.argv) > 1, "Please provide input path"
    run_folder = sys.argv[1]
    assert os.path.exists(run_folder), f"Input path must exist ({run_folder})"

    plot_occupancy(run_folder)
