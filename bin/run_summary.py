#!/usr/bin/env python3
import json
import logging
import os
import sys

import pandas as pd
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary

from column_def import ColumnDef

logger = logging.getLogger(__name__)

# Metrics summarized by read
RUN_SUMMARY_COLUMNS = [
    # Number of bases sequenced
    ColumnDef('Yield Total (G)', 'yield_g'),
    # Projected number of bases expected to be sequenced
    ColumnDef('Projected Yield (G)', 'projected_yield_g'),
    # % of passing filer clusters that aligned to the PhiX genome
    ColumnDef('% Aligned', 'percent_aligned'),
    # Calculated error rate of the reads that aligned to PhiX
    ColumnDef('Error Rate', 'error_rate'),
    # Average of the A channel intensity measured at te first cycle
    ColumnDef('Intensity C1', 'first_cycle_intensity', precision=None),
    # % of bases with a quality score of 30 or higher
    ColumnDef('% >= Q30', 'percent_gt_q30'),
    # % of clusters that can be sequenced
    ColumnDef('% Occupied', 'percent_occupied')
]

# Metrics summarized by lane
READ_SUMMARY_COLUMNS = [
    ColumnDef('Lane', 'lane', precision=None),
    # Number of tiles per lane
    ColumnDef('Tiles', 'tile_count', precision=None),
    # Density of clusters in the thousands detected by image analysis
    ColumnDef('Density', 'density', precision=None, scale=1E3),
    # Percent of clusters passing filtering
    ColumnDef('Cluster PF', 'percent_pf'),
    # The value used by RTA for the rate at which molecules in a cluster fall behind or jump ahead during a read
    ColumnDef('Legacy Phasing/Prephasing rate', ['phasing', 'prephasing'], precision=3),
    # The best-fit slope and offset of the phasing/prephasing corrections
    ColumnDef('Phasing Slope/offset', ['phasing_slope', 'phasing_offset'], precision=3),
    # See above
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
    # The calculated error rate, as determined by the PhiX alignment
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


class RunSummary:
    """
    The run summary class provides tables with basic data quality metrics summarized per lane and per read.
    """
    def __init__(self, run_folder_path):
        # Initialize interop objects
        self.run_metrics = py_interop_run_metrics.run_metrics()
        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)

        # Read from run folder
        self.run_metrics.read(run_folder_path, valid_to_load)

        # Load up summary metrics
        self.summary = py_interop_summary.run_summary()
        py_interop_summary.summarize_run_metrics(self.run_metrics, self.summary)

        # Cached result tables for subsequent calls
        self.run_summary_df = None
        self.read_summary_dfs = {}

    def get_table_html(self) -> str:
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

    def get_table_json(self) -> str:
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

    def get_run_summary(self) -> pd.DataFrame:
        """ Summarized per read """
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
        """
        Summarized per lane for a given read
        The read number starts at 0
        """
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
        """ Adds the (I) flag for indexed reads """
        read_number = self.summary.at(index).read().number()
        indexed_flag = "(I) " if self.summary.at(index).read().is_index() else ""
        return "Read {}{}".format(indexed_flag, read_number)


# If this file is being run as a standalone script
if __name__ == "__main__":

    # Get the input file path from the first argument
    assert len(sys.argv) > 1, "Please provide input path"
    input_path = sys.argv[1]
    assert os.path.exists(input_path), f"Input path must exist ({input_path})"

    # Parse the inputs
    summary = RunSummary(input_path)

    # Extract the JSON summary
    table_json = summary.get_table_json()
    with open('run_stats.json', 'wt') as handle_out:
        handle_out.write(table_json)

    table_html = summary.get_table_html()
    with open('run_stats.html', 'wt') as handle_out:
        handle_out.write(table_html)
