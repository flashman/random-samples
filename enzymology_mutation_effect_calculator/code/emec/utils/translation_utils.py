from Bio.Alphabet import Alphabet
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Data.CodonTable import CodonTable
import pandas as pd


START_SYMBOL = "M"
STOP_SYMBOL = "*"

# NCBI Standard codon table.
# See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
STANDARD_CODON_TABLE_DATA = """
    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M------**--*----M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""


def parse_codon_table(data: str, name: str = "Custom", id: int = 0) -> CodonTable:
    """
    Parse custom string-formatted codon tables as a BioPython CodonTable object.

    The string-encoded codon table data is assumed to take the form:

        AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
      Starts = ---M------**--*----M---------------M----------------------------
      Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
      Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
      Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

    Note that the protein alphabet may include non standard characters e.g. X or $.

    To access standard/pre-defined NCBI codon tables, use the Bio.Data.CodonTable methods:
    unambiguous_dna_by_name() or unambiguous_dna_by_id().

    Args:
        data: the string-encoded codon table.
        name: name of the codon table.
        id: id of the codon table.

    Returns:
        a CodonTable object.
    """
    # Read data into dataframe.
    df = pd.DataFrame()
    for line in data.splitlines():
        line = line.replace(" ", "")
        if not line:
            continue
        (row, val) = line.split("=")
        df[row] = list(val)

    # Add summary codon column
    df["Codon"] = df.Base1 + df.Base2 + df.Base3

    # Extract structured codon data.
    forward_table = {
        row.Codon: row.AAs for (i, row) in df[df.AAs != STOP_SYMBOL].iterrows()
    }
    start_codons = df[df["Starts"] == START_SYMBOL]["Codon"].tolist()
    stop_codons = df[df["Starts"] == STOP_SYMBOL]["Codon"].tolist()

    # Create custom protein alphabet to accommodate non-standard encodings.
    class CustomProteinAlphabet(Alphabet):
        letters = "".join(sorted(set(forward_table.values())))
        size = 1

    # Create the CodonTable object.
    codon_table = CodonTable(
        nucleotide_alphabet=IUPACUnambiguousDNA(),
        protein_alphabet=CustomProteinAlphabet(),
        forward_table=forward_table,
        start_codons=start_codons,
        stop_codons=stop_codons,
    )
    codon_table.id = id
    codon_table.names = [name]

    return codon_table
