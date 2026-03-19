import json
from pathlib import Path
import os
import sys
from typing import List, Optional
import warnings

from Bio import AlignIO, SeqIO
from Bio.Data.CodonTable import CodonTable, TranslationError
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

import drawbridge
from drawbridge.models.dataset import Dataset

from emec import pipeline
from emec import generalizability_pipeline
from emec.generalizability_pipeline import GeneralizabilityParams
from emec.utils import alignment_utils as alu
from emec.utils import sequence_utils as squ
from emec.utils import statsmodels_utils as stu

from zypotions.db.sql_connector import get_sqlengine

# Default input/output paths
SCRIPT_DIR = os.path.abspath(".")
INPUT_DIR = os.path.join(SCRIPT_DIR, "input_files")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output_files")

DEFAULT_SEQUENCE_FILENAME = "sequences.fasta"
DEFAULT_DNA_SEQUENCE_FILENAME = "dna_sequences.fasta"
DEFAULT_AA_SEQUENCE_FILENAME = "aa_sequences.fasta"
DEFAULT_CODON_TABLE = "Standard"

USE_REFERENCE_STRAIN = "use_reference_strain"
USE_LOGICAL_PARENT_STRAIN = "use_logical_parent_strain"

# Suppress biopython/numpy warnings. They so ugly.
warnings.simplefilter("ignore")


class MutationEffectCalculatorSession:
    """
    Session object for the MutationEffectCalculatorInterface.

    This object manages all data IO operations as well as interfacing with the core analysis
    pipeline.
    """

    def __init__(
        self, input_dir: str = INPUT_DIR, output_dir: str = OUTPUT_DIR
    ) -> None:
        """
        Initialize the session object with input and output directory paths.

        Args:
            input_dir: full path to the input directory
            output_dir: full path to the output directory
        """
        # Store data IO locations
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.dataset_id = None

        # Store analysis artifacts.
        self.alignment = None
        self.attributes_df = None
        self.results = None

        # Save db engine.
        self.engine = None

        # Set codon table.
        # See Seq.translate for valid options.
        self.codon_table = DEFAULT_CODON_TABLE

    @property
    def dataset_link(self) -> str:
        """
        Returns a link to the dataset associated with the session.
        """
        return self.lims_link("datasets", self.dataset_id)

    def load(
        self,
        performance_filename: str,
        seq_filename: str,
        reference_seq_id: int,
        feature_term: str,
        annotation_name: str,
    ) -> None:
        """
        Load sequence and performance data into the session object for subsequent analysis.

        If seq_filename is provided, load the associated sequence data.  Otherwise Use the provided
        reference_sequence_id to determine the reference protein sequence and use LIMS changes table
        payload sequence of strains indicated by the performance data file to determine the mutated
        protein sequences.

        If a feature_term and/or annotation_name is provided, queried sequence data will be
        filtered to include the sub-sequence matching these criterion.

        Args:
            performance_filename: filename of csv containing strain performance data.
            seq_filename: filename of fasta file containing sequence data.
            reference_sequence_id: sample, strian, or dna component of the reference sequence.
            feature_term: optionally filter the payload sequence by annotation feature_term.
            annotation_name: optionally filter the payload sequence by annotation name.

        Raises:
            an Exception if gap characters are introduced by the sequence alignment step.
        """
        # Load attributes
        # For now, we assume the first column contains sequence identifiers.
        self.attributes_df = self.read_csv(performance_filename, index_col=0)
        # Make sure the index is a string so that we can join with sequence data.
        self.attributes_df.index = self.attributes_df.index.astype("str")

        # Load sequence data.
        if not seq_filename and not reference_seq_id:
            raise ValueError(
                "Must provide either a reference sequence id or sequence file name."
            )

        # Query sequence data.
        if reference_seq_id and not seq_filename:
            # Name sequences file after the performance data filename.
            seq_filename_preface = Path(performance_filename).stem
            seq_filename = "_".join(
                [seq_filename_preface, DEFAULT_DNA_SEQUENCE_FILENAME]
            )
            # Query sequences but do not translate.
            self.query_sequences(
                self.attributes_df.index,
                reference_seq_id,
                feature_term,
                annotation_name,
                translate=False,
                directory=self.input_dir,
                filename=seq_filename,
            )

        # Load sequences and translate.
        seqs = self.read_sequences(
            seq_filename, translate=True, directory=self.input_dir
        )

        # Align sequences.  The first sequence is assumed to be the reference.
        # TODO(flash): Alignment parameters are currently hard-coded to avoid gap
        # characters. Possibly expose alignment params to the front end to give the user greater
        # control?
        self.alignment = alu.align_multiple_sequences_mafft(
            seqs[0], seqs[1:], op=5.0, lop=5.0
        )

        # We currently only expect missense mutations, so the alignment should be gap-less".
        # TODO(Flash): We probably want to parameterize this check?
        if any(len(self.alignment[0]) != len(s) for s in seqs):
            raise Exception("Gap characters detected in the sequence alignment.")

        # Save clustal formatted alignment to the output directory for debugging.
        fn = os.path.join(self.output_dir, "alignment.aln")
        with open(fn, "w") as f:
            AlignIO.write(self.alignment, f, "clustal")

        # Save alignment summary to the output directory for debugging.
        fn = os.path.join(self.output_dir, "alignment.csv")
        alignment_df = alu.msa_to_mutations_df(self.alignment)
        # Try to attach codon info.
        # This should fail when the input sequence data is already translated.
        try:
            dna_seqs = self.read_sequences(
                seq_filename, translate=False, directory=self.input_dir
            )
            codon_df = alu.reverse_translate_mutations_df(alignment_df, dna_seqs)
            alignment_df["target_seq_codon"] = codon_df["target_seq"]
            alignment_df["query_seq_codon"] = codon_df["query_seq"]

            # Add super special codon-based mutation name.
            # TODO(flash): This should probably go somewhere else. Shrug.
            alignment_df["name_codon"] = (
                alignment_df["target_seq_codon"]
                + alignment_df.target_position.astype(str)
                + alignment_df["query_seq_codon"]
            )

            # Save codon usage statistics for later.
            codon_usage_fn = os.path.join(self.output_dir, "codon_usage.csv")
            codon_usage_df = alu.get_codon_usage_df(
                alignment_df,
                mutation_col="name",
                codon_col=["name_codon", "query_seq_codon"],
            )
            codon_usage_df.to_csv(codon_usage_fn)

        except ValueError:
            pass

        alignment_df.to_csv(fn, index=False)

        # Save custom codon table to the output directory for dubbuging.
        # TODO(flash): For each of use, it would probably be nice to include the table using the
        # input format.
        if isinstance(self.codon_table, CodonTable):
            fn = os.path.join(self.output_dir, "codon_table.txt")
            with open(fn, "w") as f:
                f.write(str(self.codon_table))

    def run(self, perf_colname: str, pipeline_params_fn: str) -> None:
        """
        Run the analysis pipeline against performance data associated with the given column of the
        attributes file.

        Analysis results are saved to the session object and written to the output directory.

        Args:
            perf_colname: Name of performance data column in the sequence attributes dataframe.
            pipeline_params_filename: filename of json file containing pipeline parameters.
        """
        # Double check the performance file column name.
        if perf_colname not in self.attributes_df:
            valid_cols_str = "\n".join(list(self.attributes_df.columns))
            raise ValueError(
                f"Column {perf_colname} is not in the provided performance file. "
                f" Specify one of:\n {valid_cols_str}"
            )

        # Load the pipeline parameters from file.
        self.pipeline_params_fn = pipeline_params_fn
        self.pipeline_params = self.read_json(pipeline_params_fn)[
            pipeline.PIPELINE_PARAMS_KEY
        ]

        # Run the pipeline.
        self.results, self.hyperparameters, self.elastic_net = pipeline.run_pipeline(
            self.alignment, self.attributes_df[perf_colname], self.pipeline_params
        )

        # Save the fitted model to the output directory.
        fn = os.path.join(self.output_dir, "model.pickle")
        self.results.save(fn)

        # Save the fitted model predictions to the output directory.
        fn = os.path.join(self.output_dir, "model_predictions.csv")
        model_predictions_df = stu.get_model_predictions_df(self.results)
        model_predictions_df.to_csv(fn, index=False)

        # Save the model summary to the output directory.
        fn = os.path.join(self.output_dir, "model_summary.csv")
        model_summary_df = stu.get_model_summary(self.results, self.hyperparameters)
        model_summary_df.to_csv(fn)

        # Save the model coefficients to the output directory.
        fn = os.path.join(self.output_dir, "model_coefficients.csv")
        model_coefficients_df = stu.get_model_coefficients(
            self.results, drop_na=False, p_threshold=None
        )

        # Try to attach precomputed codon usage stats to model coefficient df.
        # This should fail when the input sequence data is already translated.
        try:
            # Load precomputed codon stats.
            codon_usage_fn = os.path.join(self.output_dir, "codon_usage.csv")
            codon_usage_df = pd.read_csv(codon_usage_fn, index_col="name")

            # Attach columns of interest to coefficients df.
            model_coefficients_df["codon_feature"] = codon_usage_df["name_codon_mode"]
            model_coefficients_df["codon_frequencies"] = codon_usage_df[
                "query_seq_codon_freq"
            ]
        except FileNotFoundError:
            pass

        model_coefficients_df.to_csv(fn)

    def run_generalizability_analysis(self, p_threshold: float = 0.001) -> pd.DataFrame:
        """
        Analyze how model fit and model coefficients change with the size and mutational coverage
        of the input data set. See emec.generalizability_pipeline for additional information.

        Args:
            p_threshold: filter which coefficients are plotted based on significance. Note that
            this threshold is applied to coefficient results of the main pipeline.

        Returns:
            a dataframe of model fit results for each iteration.
            a dataframe of model coefficients for each iteration.

        Raises:
            RuntimeError if this method is called before `run`.

        """

        if self.results is None:
            raise RuntimeError("Please call `self.run()` before using this method.")

        # Load the pipeline parameters from file.
        self.generalizability_params = GeneralizabilityParams(
            **self.read_json(self.pipeline_params_fn).get(
                generalizability_pipeline.GENERALIZABILITY_PARAMS_KEY, {}
            )
        )

        # Run analysis.
        fits_df, coefs_df, counts_df = generalizability_pipeline.run_pipeline(
            self.results, self.generalizability_params, self.pipeline_params
        )

        # Get coefficients displayed in main pipeline results.
        display_coefficients = stu.get_model_coefficients(
            self.results, drop_na=False, p_threshold=p_threshold
        ).index

        # Plot and save results.
        fn = self.output_dir + "/library_size_vs_model_fit.svg"
        generalizability_pipeline.plot_model_fit_vs_library_size(fits_df, saveas=fn)

        fn = self.output_dir + "/library_size_vs_feature_counts.svg"
        generalizability_pipeline.plot_feature_counts_vs_library_size(
            counts_df, saveas=fn
        )

        fn = self.output_dir + "/library_size_vs_model_hyperparameters.svg"
        generalizability_pipeline.plot_model_hyperparameters_vs_library_size(
            fits_df, saveas=fn
        )

        fn = self.output_dir + "/feature_counts_vs_model_coefficients.svg"
        generalizability_pipeline.plot_feature_counts_vs_model_coefficients(
            coefs_df, counts_df, display_coefficients, drop_zeros=True, saveas=fn
        )

        fn = self.output_dir + "/library_size_vs_model_coefficients.svg"
        generalizability_pipeline.plot_model_coefficients_vs_library_size(
            coefs_df, display_coefficients, drop_zeros=False, saveas=fn
        )

        # Return results dataframes.
        return fits_df, coefs_df, counts_df

    def save(self, use_existing_dataset: bool = True) -> None:
        """
        Save the current session (user_inputs, input_dir, output_dir) to a dataset.

        Args:
            use_existing_dataset: If true, re-use the dataset associated with the session if one
                exists. Otherwise create a new one.
        """
        self.create_dataset(use_existing=use_existing_dataset)
        self.copy_directory_to_dataset(self.input_dir, self.dataset_id)
        self.copy_directory_to_dataset(self.output_dir, self.dataset_id)

    def read_csv(
        self,
        filename: str,
        index_col: Optional[str] = None,
        directory: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Read a csv file from the given directory.

        By default, use the session's input directory.

        Args:
            filename: the name of the file to read.
            index_col: optionally, specify a column to use.
            directory: optionally, specify a directory from which to read.

        Returns:
            a DataFrame
        """
        directory = directory or self.input_dir
        fn = os.path.join(directory, filename)
        df = pd.read_csv(fn, index_col=index_col)
        return df

    def read_json(self, filename: str, directory: Optional[str] = None) -> dict:
        """
        Read a json file from the given directory.

        By default, use the session's input directory.

        Args:
            filename: the name of the file to read.
            directory: optionally, specify a directory from which to read.

        Returns:
           a dictionary.
        """
        directory = directory or self.input_dir
        fn = os.path.join(directory, filename)
        with open(fn, "r") as f:
            return json.load(f)

    def read_sequences(
        self,
        filename: str,
        fileformat: str = "fasta",
        translate: bool = True,
        directory: str = None,
    ) -> List[SeqRecord]:
        """
        Read sequences from the given directory.

        Args:
            filename: the name of the file to read.
            fileformat: the sequence file format.
            translate: if True, translate sequence data amino acid sequences.
            directory: specify a directory from which to read.

        Returns:
            a list of SeqRecord objects.
        """
        directory = directory or self.input_dir
        fn = os.path.join(directory, filename)
        seqs = [s for s in SeqIO.parse(fn, fileformat)]

        if translate:
            try:
                seqs = [
                    s.translate(
                        id=True, name=True, description=True, table=self.codon_table
                    )
                    for s in seqs
                ]
            except TranslationError:
                pass

        return seqs

    def query_sequences(
        self,
        strain_ids: List[str],
        reference_seq_id: str,
        feature_term: str,
        annotation_name: str = None,
        translate: bool = True,
        directory: Optional[str] = None,
        filename: Optional[str] = DEFAULT_SEQUENCE_FILENAME,
    ) -> List[SeqRecord]:
        """
        Query reference and mutant sequence data from LIMS.

        If reference_seq_id takes the magic value, 'use_logical_parent_strain' the reference
        sequence will be inferred from the (unique) logical parent associated with the given
        strains.

        If reference_seq_id takes the magic value, 'use_reference_strain', the reference sequence
        will be inferred from the reference strain associated with the given strains.

        Optionally save sequence data as a fasta file to the given directory.

        Args:
            strain_ids: query sequence data for the given list of strain ids
            reference_seq_id: ZId of the sample, strain, or dna component associated with the
                reference sequence. Or one of the string values: `use_reference_strain`,
               `use_logical_parent_strain`.
            feature_term: optionally filter the payload sequence by annotation feature_term.
            annotation_name: optionally filter the payload sequence by annotation name.
            translate: if True, translate sequence data amino acid sequences.
            directory: optionally, specify a directory to save the sequence data.
            filename: optionally, specify a filename to use when saving the sequence data.

        Returns:
            the queried sequence data.
        """
        # Ensure drawbridge session is set.
        # Create SQL engine using drawbridge
        drawbridge.api.sync_session()
        self.engine = get_sqlengine(
            drawbridge.DEFAULT_SESSION.instance.name, "derived_data"
        )

        # Get reference sequence.
        # First check if the reference sequence id is one of the magic values.
        if reference_seq_id == USE_REFERENCE_STRAIN:
            reference_seq_id = squ.get_changes_table_reference_strain_id(
                strain_ids, self.engine
            )
        elif reference_seq_id == USE_LOGICAL_PARENT_STRAIN:
            reference_seq_id = squ.get_changes_table_logical_parent_strain_id(
                strain_ids, self.engine
            )

        # Load the reference sequence.
        # TODO(flash): Used actual sequence name.
        ref_seq_val = squ.get_sequence(reference_seq_id, feature_term, annotation_name)
        ref_seq = SeqRecord(
            Seq(ref_seq_val), id=str(reference_seq_id), description="reference_sequence"
        )

        # Get mutant sequences.
        strain_ids = set(strain_ids) - {reference_seq_id}
        seqs = squ.get_changes_table_payload_sequences(
            strain_ids,
            self.engine,
            feature_term,
            annotation_name,
            deduplicate_strategy="last",
        )

        # Combine reference and queried sequences.
        seqs = [ref_seq] + seqs

        # Optionally translate to Amino Acid sequences.
        if translate:
            seqs = [
                s.translate(
                    id=True, name=True, description=True, table=self.codon_table
                )
                for s in seqs
            ]

        # Optionally save sequences as fasta.
        if directory:
            fn = os.path.join(directory, filename)
            with open(fn, "w") as f:
                SeqIO.write(seqs, f, "fasta")

        return seqs

    def create_dataset(self, use_existing: bool = True) -> None:
        """
        Create a new LIMS dataset and attach it to the session.

        Args:
            use_existing: if True, only create a new dataset when one isn't associated with the
                session already.

        """
        # Ensure that the drawbridge session is synced.
        drawbridge.api.sync_session()

        if not use_existing or self.dataset_id is None:
            dataset = Dataset("Enzyme Mutation Effects Calculator Results", create=True)
            self.dataset_id = dataset.id

    def copy_directory_to_dataset(self, directory: str, dataset_id: int) -> bool:
        """
        Copy contents of the given directory to the given dataset.

        Args:
            directory: path of directory to copy.
            dataset: copy data to this dataset.

        Returns:
            true on success.
        """
        dataset = Dataset.find(dataset_id)
        for fn in os.listdir(directory):
            path = os.path.join(directory, fn)
            if os.path.isdir(path):
                continue
            if fn.startswith("."):
                continue
            dataset.upload(path)

        return True

    @staticmethod
    def lims_link(model: str, model_id: int) -> str:
        """
        Generate a lims url for a particular instance of a LIMS data model for the current LIMS
        environment.

        Example::

            > lims_link("analytic_summaries", 123)
            https://lims.zymergen.net/analytic_summaries/123

        Arguments:
            model: the data model name, .
            model_id: Id of model instance.

        Returns:
            a LIMS url to the instance.

        Raises:
            an Exception if a drawbridge exception has not been established.
        """
        if not drawbridge.DEFAULT_SESSION:
            raise Exception("Drawbridge not connected to LIMS instance.")

        # Hack the zwork base url to get the LIMS base url.
        # TODO(flash): Update drawbridge to create these urls instead.
        zwork_url = drawbridge.DEFAULT_SESSION.instance.zwork_api_url
        lims_url = zwork_url.replace("-jobs", "")

        return f"{lims_url}/{model}/{model_id}"
