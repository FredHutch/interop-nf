import json
import unittest

from bin.run_summary import RunSummary
from tests.test_helper import miseq_demo_path


class RunSummaryTest(unittest.TestCase):
    """ See readme in tests/data for downloading test dataset """

    def test_parse_miseq_data_json(self):
        expected_reads = 2
        expected_lanes = 1
        summary = RunSummary(miseq_demo_path)

        table_json = summary.get_table_json()

        parsed = json.loads(table_json)
        run_summary = parsed['runSummary']
        reads = parsed['reads']

        self.assertIsNotNone(parsed)

        # Check run summary
        self.assertEqual(expected_reads + 2, len(run_summary))
        for level in run_summary:
            self.assertIsNotNone(level['Level'])
            self.assertTrue(level['Yield Total (G)'] > 0)
            self.assertTrue(level['% Aligned'] > 0)

        # Check reads section, length matches number of reads and number of lanes.
        # Yield is present
        self.assertEqual(expected_reads, len(reads))
        for read in reads:
            self.assertEqual(expected_lanes, len(read))
            self.assertIsNotNone(read[0]['Yield'])
            self.assertEqual(1, read[0]['Lane'])

    def test_parse_miseq_data_html(self):
        # Run summary plus two read tables
        expected_number_of_tables = 3
        summary = RunSummary(miseq_demo_path)

        table_html = summary.get_table_html()

        self.assertIsNotNone(table_html)
        self.assertTrue('Run Quality Summary' in table_html)

        number_of_tables = table_html.count('<table')
        self.assertEqual(expected_number_of_tables, number_of_tables)
