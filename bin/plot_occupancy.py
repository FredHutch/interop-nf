#!/usr/bin/env python3

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from interop import py_interop_run_metrics, py_interop_run, py_interop_table


def plot_occupancy(run_folder: str, output_jpg_prefix="occupancy"):
    """
    To optimize loading concentrations on the NovaSeq platform, the % Occupied and % Pass Filter
    metrics can be plotted to determine if a run was underloaded, optimally loaded, or overloaded.

    More information:
    https://support.illumina.com/bulletins/2020/03/plotting---occupied-by---pass-filter-to-optimize-loading-concent.html
    """

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
                [f"{column.name()} ({subname})" for subname in column.subcolumns()])
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

    # Skip if there is no data (% Occupied only available on NovaSeq)
    if df.shape[0] == 0 or "% Occupied" not in df:
        # Stop
        print("Occupancy plot skipped, no data available")
        return

    x = "% Occupied"
    y = "% Pass Filter"
    hues = ["Tile", "Lane", "Cycle"]

    # Make a few different types of plots
    for hue in hues:
        sns.scatterplot(
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
        plt.savefig(f"{output_jpg_prefix}_{hue.lower()}.jpg", dpi=600)
        plt.close()


# If this file is being run as a standalone script
if __name__ == "__main__":

    # Get the input file path from the first argument
    assert len(sys.argv) > 1, "Please provide input path"
    input_path = sys.argv[1]
    assert os.path.exists(input_path), f"Input path must exist ({input_path})"

    plot_occupancy(input_path)
