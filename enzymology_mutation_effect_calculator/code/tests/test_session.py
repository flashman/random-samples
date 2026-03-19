import os
import shutil
import tempfile
import unittest

from emec.utils import translation_utils
from session import MutationEffectCalculatorSession
from tests.fixtures.random_proteins_simulated_small import (
    DNA_SEQUENCE_PATH,
    SEQUENCE_PATH,
    ATTRIBUTES_PATH,
    PARAMS_PATH,
)


class TestSession(unittest.TestCase):
    def setUp(self):
        # Create temp directories
        input_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, input_dir)

        output_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, output_dir)

        # Add fixtures to temp directory
        shutil.copy(DNA_SEQUENCE_PATH, input_dir)
        shutil.copy(SEQUENCE_PATH, input_dir)
        shutil.copy(ATTRIBUTES_PATH, input_dir)
        shutil.copy(PARAMS_PATH, input_dir)

        self.session = MutationEffectCalculatorSession(input_dir, output_dir)

    def test_init(self):
        self.assertTrue(self.session.input_dir)
        self.assertTrue(self.session.output_dir)

    def test_load(self):
        self.session.load("attributes.csv", "sequences.fa", None, None, None)

        self.assertLess(0, len(self.session.attributes_df))
        self.assertLess(0, len(self.session.alignment))

        assert os.path.exists(os.path.join(self.session.output_dir, "alignment.aln"))
        assert os.path.exists(os.path.join(self.session.output_dir, "alignment.csv"))

    def test_load_dna(self):
        self.session.load("attributes.csv", "dna_sequences.fa", None, None, None)

        self.assertLess(0, len(self.session.attributes_df))
        self.assertLess(0, len(self.session.alignment))

        assert os.path.exists(os.path.join(self.session.output_dir, "alignment.aln"))
        assert os.path.exists(os.path.join(self.session.output_dir, "alignment.csv"))

    def test_load_dna_custom_codon_table(self):
        table = translation_utils.parse_codon_table(
            translation_utils.STANDARD_CODON_TABLE_DATA
        )
        self.session.codon_table = table

        self.session.load("attributes.csv", "dna_sequences.fa", None, None, None)

        self.assertLess(0, len(self.session.attributes_df))
        self.assertLess(0, len(self.session.alignment))

        assert os.path.exists(os.path.join(self.session.output_dir, "alignment.aln"))
        assert os.path.exists(os.path.join(self.session.output_dir, "alignment.csv"))

    def test_run(self):
        self.session.load("attributes.csv", "sequences.fa", None, None, None)
        self.session.run("fiop", "pipeline_parameters.json")
        self.assertTrue(self.session.results)
        self.assertTrue(self.session.hyperparameters)

    def test_run_dna(self):
        self.session.load("attributes.csv", "dna_sequences.fa", None, None, None)
        self.session.run("fiop", "pipeline_parameters.json")
        self.assertTrue(self.session.results)
        self.assertTrue(self.session.hyperparameters)

    def test_read_sequences(self):
        seqs = self.session.read_sequences(
            "dna_sequences.fa", translate=False, directory=self.session.input_dir
        )
        self.assertEqual(2001, len(seqs))
        self.assertEqual(30, len(seqs[0].seq))

        seqs = self.session.read_sequences(
            "dna_sequences.fa", translate=True, directory=self.session.input_dir
        )
        self.assertEqual(2001, len(seqs))
        self.assertEqual(10, len(seqs[0].seq))

        seqs = self.session.read_sequences(
            "sequences.fa", translate=True, directory=self.session.input_dir
        )
        self.assertEqual(2001, len(seqs))
        self.assertEqual(10, len(seqs[0].seq))
