import os.path
import pathlib

FIXTURES_PATH = pathlib.Path(os.path.dirname(__file__))

DNA_SEQUENCE_PATH = FIXTURES_PATH / "dna_sequences.fa"
SEQUENCE_PATH = FIXTURES_PATH / "sequences.fa"
ATTRIBUTES_PATH = FIXTURES_PATH / "attributes.csv"
PARAMS_PATH = FIXTURES_PATH / "pipeline_parameters.json"
