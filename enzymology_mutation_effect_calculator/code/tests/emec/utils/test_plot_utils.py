import pickle
import unittest

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import emec.utils.plot_utils as pu
import emec.utils.statsmodels_utils as su

from tests.fixtures.fern_gar0722_ssl_small_2021_08_04 import MODEL_PATH


class TestPlotUtils(unittest.TestCase):
    def test_plot_library_coverage(self):
        """
        Minimally exercise the code.
        """
        msa = MultipleSeqAlignment(
            [
                SeqRecord(seq=Seq("actga"), id="0"),
                SeqRecord(seq=Seq("accga"), id="1"),
                SeqRecord(seq=Seq("aacta"), id="2"),
            ]
        )

        f = pu.plot_library_coverage(msa)

        self.assertTrue(f)

    def test_plot_mutation_co_occurence(self):
        """
        Minimally exercise the code.
        """
        msa = MultipleSeqAlignment(
            [
                SeqRecord(seq=Seq("actga"), id="0"),
                SeqRecord(seq=Seq("accga"), id="1"),
                SeqRecord(seq=Seq("aacta"), id="2"),
            ]
        )

        f = pu.plot_mutation_co_occurence(msa)

        self.assertTrue(f)

    def test_plot_sequence_distance(self):
        """
        Minimally exercise the code.
        """
        msa = MultipleSeqAlignment(
            [
                SeqRecord(seq=Seq("actga"), id="0"),
                SeqRecord(seq=Seq("accga"), id="1"),
                SeqRecord(seq=Seq("aacta"), id="2"),
            ]
        )

        f = pu.plot_sequence_distance(msa)

        self.assertTrue(f)

    def test_plot_model_coefficients(self):
        """
        Minimally exercise the code.
        """
        with open(MODEL_PATH, "rb") as f:
            results = pickle.load(f)

        f = pu.plot_model_coefficients(results)

        self.assertTrue(f)

    def test_plot_model_coefficient_interactions(self):
        """
        Minimally exercise the code.
        """
        with open(MODEL_PATH, "rb") as f:
            results = pickle.load(f)

        f = pu.plot_model_coefficient_interactions(results)

        self.assertTrue(f)

    def test_plot_model_fit(self):
        """
        Minimally exercise the code.
        """
        with open(MODEL_PATH, "rb") as f:
            results = pickle.load(f)

        f = pu.plot_model_fit(results)

        self.assertTrue(f)

    def test_plot_model_coefficient_support(self):
        """
        Minimally exercise the code.
        """
        with open(MODEL_PATH, "rb") as f:
            results = pickle.load(f)

        f = pu.plot_model_coefficient_support(results)

        self.assertTrue(f)

    def test_get_colors(self):
        """
        Minimally exercise the code.
        """
        with open(MODEL_PATH, "rb") as f:
            results = pickle.load(f)

        df = su.get_model_coefficients(results)

        c = pu._get_colors(df.coef)
        self.assertEqual(len(c), len(df))

        c = pu._get_colors(df.coef, palette=pu.DEFAULT_COEF_PALETTE)
        self.assertEqual(len(c), len(df))

        c = pu._get_colors(df.coef, keys=df.index)
        self.assertEqual(list(c), list(df.index))

        c = pu._get_colors(df.coef, center=0)
        self.assertEqual(len(c), len(df))
