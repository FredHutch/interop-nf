import unittest
from pathlib import Path

import matplotlib.pyplot as plt

from bin.plot_occupancy import plot_occupancy
from bin.plot_percent_base import plot_percent_base
from bin.plot_tile_intensity import plot_tile_intensity
from tests.test_helper import miseq_demo_path


class PlotTests(unittest.TestCase):

    def test_plot_percent_base(self):
        plt.figure()
        plot_percent_base(miseq_demo_path)
        self.assertTrue(Path('percent_base.svg').exists())

    def test_plot_tile_intensity(self):
        plt.figure()
        plot_tile_intensity(miseq_demo_path)
        self.assertTrue(Path('max_intensity_1.svg').exists())

    def test_plot_occupancy_should_not_created_on_demo(self):
        """ MiSeq shouldn't have this data """
        plt.figure()
        plot_occupancy(miseq_demo_path)
        self.assertFalse(Path('occupancy_tile.jpg').exists())
        self.assertFalse(Path('occupancy_cycle.jpg').exists())
        self.assertFalse(Path('occupancy_read.jpg').exists())
