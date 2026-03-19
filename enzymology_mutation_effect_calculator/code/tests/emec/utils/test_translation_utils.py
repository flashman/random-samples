from Bio.Data import CodonTable
from Bio.Seq import Seq
import unittest

from emec.utils import translation_utils


# Modified versoin of the NCBI Standard codon table:
#   CTG -> O
#   GTG -> z
CUSTOM_CODON_TABLE_DATA = """
    AAs  = FFLLSSSSYY**CC*WLLLOPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVZAAAADDEEGGGG
  Starts = ---M------**--*----M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""


# Modified versoin of the NCBI Standard codon table using non-alphabetic characters:
#   CTG -> @
#   GTG -> &
SYMBOLS_CODON_TABLE_DATA = """
    AAs  = FFLLSSSSYY**CC*WLLL@PPPPHHQQRRRRIIIMTTTTNNKKSSRRVVV&AAAADDEEGGGG
  Starts = ---M------**--*----M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""


class TestTranslationUtils(unittest.TestCase):
    def test_parse_standard(self):
        table = translation_utils.parse_codon_table(
            translation_utils.STANDARD_CODON_TABLE_DATA, "Standard, SGC0", 1
        )
        standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.assertEqual(str(standard_table), str(table))

    def test_parse_custom(self):
        table = translation_utils.parse_codon_table(CUSTOM_CODON_TABLE_DATA)

        self.assertTrue("O" in table.protein_alphabet.letters)
        self.assertTrue("Z" in table.protein_alphabet.letters)

        self.assertEqual("O", table.forward_table["CTG"])
        self.assertEqual("Z", table.forward_table["GTG"])

    def test_parse_custom_symbols(self):
        table = translation_utils.parse_codon_table(SYMBOLS_CODON_TABLE_DATA)

        self.assertTrue("@" in table.protein_alphabet.letters)
        self.assertTrue("&" in table.protein_alphabet.letters)

        self.assertEqual("@", table.forward_table["CTG"])
        self.assertEqual("&", table.forward_table["GTG"])

    def test_translate_custom(self):
        table = translation_utils.parse_codon_table(CUSTOM_CODON_TABLE_DATA)
        seq = Seq("ATGCTGGTG")
        aa = seq.translate(table)

        self.assertEqual("MOZ", str(aa))

    def test_translate_custom_symbols(self):
        table = translation_utils.parse_codon_table(SYMBOLS_CODON_TABLE_DATA)
        seq = Seq("ATGCTGGTG")
        aa = seq.translate(table)

        self.assertEqual("M@&", str(aa))
