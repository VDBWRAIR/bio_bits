from . import unittest

from bio_bits import plot_muts

class TestExtractDate(unittest.TestCase):
    def test_valid_header(self):
        h = '>foo____2000_01_01'
        r = plot_muts.extract_date(h)
        self.assertEqual(r.year, 2000)
        self.assertEqual(r.month, 1)
        self.assertEqual(r.day, 1)

    def test_missing_underscores(self):
        h = '>foo'
        self.assertRaises(
            plot_muts.InvalidFastaIdentifier, plot_muts.extract_date, h
        )
