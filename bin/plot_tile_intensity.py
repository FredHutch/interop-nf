#!/usr/bin/env python3

from interop import py_interop_run_metrics, py_interop_run
import os
import pandas as pd
import seaborn as sns
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_tile_intensity(run_folder):

    # Initialize interop objects
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    valid_to_load[py_interop_run.Extraction] = 1

    # Read from the run folder
    run_metrics.read(run_folder, valid_to_load)

    # Get the extraction metrics
    extraction_metrics = run_metrics.extraction_metric_set()
    extraction_metrics.rebuild_index(True)

    # Format a DataFrame
    df = []
    for lane in extraction_metrics.lanes():

        print(f"Processing lane {lane}")

        for tile in extraction_metrics.tile_numbers_for_lane(lane):

            if not isinstance(tile, int):
                continue
            
            print(f"Processing tile {tile}")

            for cycle in range(extraction_metrics.max_cycle()):

                print(f"Processing cycle {cycle}")

                try:

                    extraction_metric = extraction_metrics.get_metric(
                        lane,
                        tile,
                        cycle
                    )

                except:

                    continue

                for channel in range(extraction_metrics.channel_count()):

                    df.append(dict(
                        channel=channel,
                        lane=lane,
                        tile=tile,
                        cycle=cycle,
                        max_intensity=extraction_metric.max_intensity(channel)
                    ))
        
    df = pd.DataFrame(df)

    # Make the plot
    with PdfPages("max_intensity.pdf") as pdf:

        # Iterate over lanes
        for lane, lane_df in df.groupby('lane'):

            # Plot the change in max intensity over cycles
            g = sns.lineplot(
                data=lane_df,
                x='cycle',
                y='max_intensity',
            )

            # Set a title
            plt.title(f"Lane {lane}")

            # Save to the PDF
            pdf.savefig()
            plt.close()


# If this file is being run as a standalone script
if __name__ == "__main__":

    # Get the input file path from the first argument
    assert len(sys.argv) > 1, "Please provide input path"
    run_folder = sys.argv[1]
    assert os.path.exists(run_folder), f"Input path must exist ({run_folder})"

    plot_tile_intensity(run_folder)
