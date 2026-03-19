from collections import defaultdict
import os
import random

from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
import pandas as pd

from .alignment_utils import (
    align_multiple_sequences_mafft,
    get_mutations,
    msa_to_mutations_df,
)


class RandomMutagenesisSimulation:
    """
    Class to simulation random mutagenesis experiments.

    Example::

        # Configure the simulation such that mutants contain 3 mutations
        # and such that there are 5 distinct beneficial mutations.
        mutation_params = {"n_mutations": 3}
        fiop_params = {"n_beneficial": 5, "benefit_mu": 2.0, "benefit_sigma": 0.05}

        # Run the simulation.
        sim = RandomMutagenesisSimulation(10, 2000, mutation_params, fiop_params)
        sim.simulate_sequences()
        sim.simulate_fiops_v1()

        # Save the simulation to a directory.
        sim.save('simulation/')

    TODO(flash): Create abstract base class to support other types of simulations.
    """

    def __init__(
        self,
        seq_length: int,
        library_size: int,
        mutation_params: dict = None,
        fiop_params: dict = None,
        reference_sequence: SeqRecord = None,
        alphabet: IUPAC.Alphabet = IUPAC.protein,
    ):
        """
        Initialize the random mutagenesis simulator.

        Args:
            seq_length: length of simulated sequences.
            library_size: the number of mutated sequences
            mutation_params: dictionary of parameters controlling mutant sequence generation.
                These parameters are specifically passed to the `mutate_sequence` method.
            fiop_params: dictionary of parameters controlling how fiop values are calculated.
                These parameters are specifically passed to the `simulate_fiop_v1` method.
            reference_sequence: optional sequence to use as the reference. If none is specified, a
                random sequence will be created.
            alphabet: specify a specific alphabet to use in the simulation.
        """

        # Simulation parameters.
        self.alphabet = alphabet
        self.seq_length = seq_length
        self.library_size = library_size
        self.mutation_params = mutation_params or defaultdict(lambda: None)
        self.fiop_params = fiop_params or defaultdict(lambda: None)

        # Simulation storage.
        self.reference_sequence = reference_sequence
        self.mutant_sequences = None
        self.alignment = None
        self.attributes_df = None

    def simulate_sequences(self) -> None:
        """
        Simulate random mutants of a random reference sequence using the `mutate_sequence` method
        and align them to the reference.
        """
        # Init ref sequence data.
        if self.reference_sequence is None:
            self.reference_sequence = random_sequence(
                self.seq_length, id="-1", alphabet=self.alphabet
            )

        # Init mutant sequence data.
        self.mutant_sequences = [
            mutate_sequence(self.reference_sequence, id=str(i), **self.mutation_params)
            for i in range(self.library_size)
        ]

        # Create alignment from sequences.
        # Alignment params are hard coded to avoid gap characters.
        self.alignment = align_multiple_sequences_mafft(
            self.reference_sequence, self.mutant_sequences, op=5.0, lop=5.0
        )

    def simulate_fiops_v1(self) -> None:
        """
        Simulate fiop data using the `simulate_fiop_v1` method.

        See `simulate_fiop_v1` documentation for additional information.
        """
        if self.reference_sequence is None or self.mutant_sequences is None:
            raise Exception("Must run 'simulate_sequences' before calling this method.")

        self.attributes_df, self.model_df = simulate_fiop_v1(
            self.reference_sequence, self.mutant_sequences, **self.fiop_params
        )

    def save(self, directory: str) -> None:
        """
        Save the simulated to the given directory.

        The saved simulation includes:
        1. the aligned sequences, as a fasta file
        2. simulated performance data, as a csv
        3. model parameters, as a csv.
        """
        fn = os.path.join(directory, "sequences.fa")
        with open(fn, "w") as f:
            AlignIO.write(self.alignment, f, "fasta")

        fn = os.path.join(directory, "attributes.csv")
        self.attributes_df.to_csv(fn)

        fn = os.path.join(directory, "model.csv")
        self.model_df.to_csv(fn)


def random_sequence(n, alphabet=IUPAC.unambiguous_dna, **kwargs):
    """
    Create a random sequence using the given alphabet.

    Args:
        n: the length of the random sequence.
        alphabet: the alphabet from which the random sequence is drawn.
        **kwargs: any additional arguments to forward to the SeqRecord constructor.

    Returns:
        a SeqRecord object.
    """
    seq = Seq("".join(random.choices(list(alphabet.letters), k=n)), alphabet=alphabet)
    return SeqRecord(seq, **kwargs)


def mutate_sequence(sequence, n_mutations, alphabet=None, positions=None, **kwargs):
    """
    Mutate a given sequence by swapping characters at exactly n positions.

    If alphabet is specified use it to source mutations, otherwise use the alphabet of the provided
    parent sequence.

    If positions are provided, only perform mutations at those positions.

    Args:
        sequence: a SeqRecord object for the sequence being mutated.
        n_mutations: the number of positions in the sequence to mutate.
        alphabet: optional Bio.Alphabet instance to use when creating mutations.  If not specified,
            use the alphabet of the given sequence.
        positions: optional positions to target for mutation.  If not specified, use any position.
        library associated with seq.
        **kwargs: optional arguments passed to the newly created SeqRecord object.

    Returns:
        a SeqRecord object with mutated sequence..
    """
    s = sequence.seq
    alphabet = alphabet or s.alphabet

    # Shuffle positions and take first n.
    positions = positions or list(range(len(s)))
    random.shuffle(positions)
    mutate_positions = set(positions[:n_mutations])

    # Mutate randomly selected positions.
    m = Seq(
        "".join(
            random.choice(list(set(alphabet.letters) - {str(a)}))
            if i in mutate_positions
            else str(a)
            for i, a in enumerate(str(s))
        ),
        alphabet=alphabet,
    )

    return SeqRecord(m, **kwargs)


def simulate_fiop_v1(
    reference,
    sequences,
    n_beneficial=10,
    mu=1.0,
    sigma=0.05,
    benefit_mu=1.7,
    benefit_sigma=0.5,
):
    """
    Simulate FIOP data for a collection of mutated sequences.

    Simulated data is generated according to the following procedure:
        1. Sample FIOP values for each mutant sequence from N(1, sigma)
        2. Select n beneficial mutations
        3. Assign a benefit to each of these mutations.  Benefits are sampled from
           N(benefit_mu, benefit_sigma)
        4. For each sequence containing beneficial mutations, apply all benefits to the base FIOP
           value by multiplying them together.

    Args:
        reference: a SeqRecord object corresponding to the base sequence.
        sequences: a list of SeqRecord objects corresponding to the mutated sequences.
        n_beneficial: the number of beneficial mutations in the simulation.
        mu: the mean of sequences lacking beneficial_mutations.
        sigma: the variance of sequences lacking beneficial mutations.
        benefit_mu: mean benefit of beneficial mutations.
        benefit_sigma: variance of benefit of beneficial mutations.

    Returns:
        df: a dataframe of simulated piop value, indexed on sequence id.
        model_df: a dataframe summarizing the beneficial_mutations.
    """

    # Identify all unique mutations relative to the reference.
    ali = align_multiple_sequences_mafft(reference, sequences)
    mut_df = msa_to_mutations_df(ali)

    # Assign normal performance to all sequences
    df = pd.DataFrame(
        np.random.normal(mu, sigma, len(ali)),
        index=[r.id for r in ali],
        columns=["fiop"],
    )

    # Fiop should never be less than 0.
    df["fiop"] = np.maximum(df["fiop"], 0)

    # Use mutant sequence to generate beneficial mutations and associated benefits.
    beneficial_mutations = get_mutations(
        reference, mutate_sequence(reference, n_beneficial, id="mutant")
    )
    effects = np.maximum(0, np.random.normal(benefit_mu, benefit_sigma, n_beneficial))

    # Apply effect to effected sequences.
    for bm, eff in zip(beneficial_mutations, effects):
        improved = mut_df[mut_df.name == bm.name].query_id.unique()
        df.loc[improved, "fiop"] = df.loc[improved, "fiop"] * eff

    # Prep beneficial mutations for output.
    model_df = pd.DataFrame(beneficial_mutations)
    model_df["effect"] = effects

    return df, model_df
