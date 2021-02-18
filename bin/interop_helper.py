#!/usr/bin/env python3

import json
import logging
import os
from pathlib import Path
import sys
from typing import Union, List

import math
import matplotlib.pyplot as plt
import pandas as pd
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot
from interop.py_interop_summary import run_summary
from matplotlib.figure import Figure

logger = logging.getLogger(__name__)


class ColumnDef(object):
    """ Represents a data value from the Illumina InterOp output """
    def __init__(self, name: str, field: Union[str, List[str]], **kwargs):
        self.name = name
        if type(field) is str:
            self.fields = [field]
        else:
            self.fields = field
        self.precision = kwargs.get('precision', 2)
        self.scale = kwargs.get('scale', 1)

    def get_value_from_row(self, row) -> str:
        """
        The attributes are a proxy of the c implementation, the values are actually a method getter
        This method works with values that have mean/std deviations, ranges, and multiple values
        """

        if len(self.fields) > 1:
            # Special case where there are multiple fields
            values = []
            for field in self.fields:
                values.append(str(self._format(getattr(row, field)().mean())))
            return ' / '.join(values)

        value_method = getattr(row, self.fields[0])
        if value_method is None:
            return 'N/A'

        value = value_method()

        if hasattr(value, 'mean'):
            # Special case when the value has a mean, add standard deviation
            return f'{self._format(value.mean())} +/- {self._format(value.stddev())}'
        elif hasattr(value, 'error_cycle_range'):
            # Special case when there is a cycle range
            error_range = value.error_cycle_range()
            first = error_range.first_cycle()
            last = error_range.last_cycle()
            if first == last:
                return self._format(first)
            return f'{self._format(first)} - {self._format(last)}'
        else:
            return self._format(value)

    def _format(self, value):
        value = value / self.scale

        if not math.isnan(value):
            value = round(value, self.precision)
        return value


RUN_SUMMARY_COLUMNS = [
    # Number of bases sequenced
    ColumnDef('Yield Total (G)', 'yield_g'),
    # Projected number of bases expected to be sequenced
    ColumnDef('Projected Yield (G)', 'projected_yield_g'),
    # Aligned to te PhiX genome
    ColumnDef('% Aligned', 'percent_aligned'),
    # Calculated error rate of the reads that aligned to PihX
    ColumnDef('Error Rate', 'error_rate'),
    # Average of the A channel intensity measured at te first cycle
    ColumnDef('Intensity C1', 'first_cycle_intensity', precision=None),
    # % of bases with a quality score of 30 or higher
    ColumnDef('% >= Q30', 'percent_gt_q30'),
    # % of clusters that can be sequenced
    ColumnDef('% Occupied', 'percent_occupied')
]

READ_SUMMARY_COLUMNS = [
    ColumnDef('Lane', 'lane', precision=None),
    ColumnDef('Tiles', 'tile_count', precision=None),
    # Density of clusters in the thousands detected by image analysis
    ColumnDef('Density', 'density', precision=None, scale=1E3),
    # Percent of clusters passing filtering
    ColumnDef('Cluster PF', 'percent_pf'),
    ColumnDef('Legacy Phasing/Prephasing rate', ['phasing', 'prephasing'], precision=3),
    ColumnDef('Phasing Slope/offset', ['phasing_slope', 'phasing_offset'], precision=3),
    ColumnDef('Prephasing Slope/offset', ['prephasing_slope', 'prephasing_offset'], precision=3),
    # Number of clusters (in millions)
    ColumnDef('Reads', 'reads', scale=1E6),
    # Number of clusters passing filtering
    ColumnDef('Reads PF', 'reads_pf', scale=1E6),
    # Number of bases with a quality score of 30 or higher
    ColumnDef('% >= Q30', 'percent_gt_q30'),
    # Yield total
    ColumnDef('Yield', 'yield_g'),
    ColumnDef('Cycles Error', 'cycle_state', precision=None),
    ColumnDef('Aligned', 'percent_aligned'),
    ColumnDef('Error', 'error_rate'),
    ColumnDef('Error (35)', 'error_rate_35'),
    ColumnDef('Error (50)', 'error_rate_50'),
    ColumnDef('Error (100)', 'error_rate_100'),
    ColumnDef('Intensity C1', 'first_cycle_intensity', precision=None),
    ColumnDef('% Occupied', 'percent_occupied'),
]

SAV_TABLE_STYLE = """
<style>
.sav-table {
  border-collapse: collapse;
}

.sav-table th, .sav-table td {
  padding: 4px;
  border: 1px solid #ddd;
  text-align: center;
}
</style>
"""


class InterOpHelper:
    """
    Wrapper for the Illumina Interop Library
    http://illumina.github.io/interop/index.html
    """
    def __init__(self, run_folder_path: str):
        self.run_folder = run_folder_path
        self.run_metrics = None
        self.run_summary_df = None
        self.read_summary_dfs = {}
        self.summary = self._load_metrics()

    def get_run_stats_table_html(self) -> str:
        logger.info('Generating run stats html')
        read_tables = ''.join([
            f"""
            <br>
            <strong>Read {read_num + 1}</strong> 
            {self.get_read_summary(read_num).to_html(index=False, border=0, classes=['sav-table'])}
            """
            for read_num in range(0, self.summary.size())
        ])
        return f"""
        {SAV_TABLE_STYLE}
        <strong>Run Quality Summary</strong>
        {self.get_run_summary().to_html(index=False, border=0, classes=['sav-table'])}
        {read_tables}
        """

    def get_run_stats_table_json(self):
        logger.info('Generating run stats json')
        summary_data = {
            'runSummary': self.get_run_summary().to_dict('records'),
            'reads': [
                self.get_read_summary(read_num).to_dict('records')
                for read_num in range(0, self.summary.size())
            ]
        }
        json_output = json.dumps(summary_data, indent=2)
        # Look into using third party JSON library if need more formatting options
        json_output = json_output.replace('NaN', 'null')
        return json_output

    def save_percent_base_plot(self, directories: List[str]) -> None:
        """
        Saves the % Base plot under the run and transfer directory as a png
        """
        base_percent_plot = self.create_percent_base_plot()

        for folder in directories:
            plot_path = os.path.join(folder, 'percent_base_plot.png')
            if Path(plot_path).exists():
                logger.info(f'{plot_path} already exists, skipping save')
                continue

            base_percent_plot.savefig(plot_path)
            logger.info(f'Saved plot to {plot_path}')

    def create_percent_base_plot(self) -> Figure:
        """
        This plots the base % across each cycle. Each line represents a different base
        """
        logger.info('Generating % base plot')
        plot_data = py_interop_plot.candle_stick_plot_data()
        options = py_interop_plot.filter_options(self.run_metrics.run_info().flowcell().naming_method())
        py_interop_plot.plot_by_cycle(self.run_metrics, "BasePercent", options, plot_data)

        figure = plt.figure()

        # Plot each base
        for base_index in range(plot_data.size()):
            line_data = plot_data.at(base_index)
            x = [line_data.at(i).x() for i in range(line_data.size())]
            y = [line_data.at(i).y() for i in range(line_data.size())]
            plt.plot(x, y, color=line_data.color(), linewidth=0.5, label=line_data.title())

        # Plot reference lines for reads
        read_vector = self.run_metrics.run_info().reads()
        for read_index in range(read_vector.size()):
            read_name = f'R{read_vector[read_index].number()}'
            cycle_start = read_vector[read_index].first_cycle()
            plt.axvline(x=cycle_start, color='purple', linestyle='--', linewidth=0.35)
            plt.text(cycle_start, plt.gca().get_ylim()[1], read_name, fontsize=8, color='purple')

        axes_data = plot_data.xyaxes()
        plt.xlabel(axes_data.x().label(), fontsize=10)
        plt.ylabel(axes_data.y().label(), fontsize=10)
        plt.title(plot_data.title(), fontsize=10)
        plt.legend()
        plt.ylim([axes_data.y().min(), axes_data.y().max()])
        plt.xlim([axes_data.x().min(), axes_data.x().max()])
        return figure

    def get_run_summary(self) -> pd.DataFrame:
        """
        The Summary table contains basic data quality metrics summarized per lane
        and per read
        """
        if self.run_summary_df is not None:
            return self.run_summary_df

        rows = [(self._get_read_display_name(read_num), self.summary.at(read_num).summary())
                for read_num in range(0, self.summary.size())]
        rows += [('Non-Indexed Total', self.summary.nonindex_summary()), ('Total', self.summary.total_summary())]

        data = {}
        for column in RUN_SUMMARY_COLUMNS:
            data[column.name] = pd.Series([column.get_value_from_row(row[1]) for row in rows],
                                          index=[row[0] for row in rows])
        df = pd.DataFrame.from_dict(data)
        df.index.name = 'Level'
        df.reset_index(inplace=True)
        self.run_summary_df = df
        return self.run_summary_df

    def get_read_summary(self, read_num: int) -> pd.DataFrame:
        """ This table repeats per lane, starts at 0 """
        if read_num in self.read_summary_dfs:
            return self.read_summary_dfs[read_num]

        if read_num < 0:
            raise RuntimeError("Read number must be greater or equal to 0")
        if read_num > self.summary.size():
            raise RuntimeError("Read number is greater than available reads")
        rows = [self.summary.at(read_num).at(lane) for lane in range(0, self.summary.lane_count())]

        data = {}
        for column in READ_SUMMARY_COLUMNS:
            data[column.name] = pd.Series([column.get_value_from_row(row) for row in rows])
        self.read_summary_dfs[read_num] = pd.DataFrame.from_dict(data)
        return self.read_summary_dfs[read_num]

    def _get_read_display_name(self, index: int) -> str:
        read_number = self.summary.at(index).read().number()
        indexed_flag = "(I) " if self.summary.at(index).read().is_index() else ""
        return "Read {}{}".format(indexed_flag, read_number)

    def _load_metrics(self) -> run_summary:
        """ Initializes the run metrics """
        self.run_metrics = py_interop_run_metrics.run_metrics()
        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
        self.run_metrics.read(self.run_folder, valid_to_load)
        summary = py_interop_summary.run_summary()
        py_interop_summary.summarize_run_metrics(self.run_metrics, summary)
        return summary


# If this file is being run as a standalone script
if __name__ == "__main__":

    # Get the input file path from the first argument
    assert len(sys.argv) > 1, "Please provide input path"
    input_fp = sys.argv[1]
    assert os.path.exists(input_fp), f"Input path must exist ({input_fp})"

    # Parse the inputs
    helper = InterOpHelper(input_fp)

    # Extract the JSON summary
    table_json = helper.get_run_stats_table_json()
    with open('qc_metrics.json', 'wt') as handle_out:
        handle_out.write(table_json)

    # Write out the plot
    helper.create_percent_base_plot().savefig("percent_base_plot.png")

    # Write out an HTML summary
    with open('run_stats.html', 'wt') as handle_out:
        handle_out.write(helper.get_run_stats_table_html())
