from collections import namedtuple
import io
import re
import tempfile
from typing import Dict, List, Union

from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import scipy.spatial.distance as dist


Mutation = namedtuple(
    "Mutation",
    field_names=[
        "target_id",
        "query_id",
        "alignment_position",
        "target_position",
        "query_position",
        "target_seq",
        "query_seq",
        "name",
        "offset",
    ],
)

DEFAULT_OFFSET = 1
GAP_SYM = "-"
COMP_COL_SORT_RE = re.compile(r"(?P<pos>\d+)")


def align_multiple_sequences_mafft(
    target_seq: SeqRecord, query_seqs: List[SeqRecord], **kwargs
) -> MultipleSeqAlignment:
    """
    Perform multiple sequence alignment of query sequences against a target/reference sequence using
    MAFFT.

    See https://biopython.org/DIST/docs/api/Bio.Align.Applications._Mafft-pysrc.html for a list of
    valid kwargs.

    Args:
        target_seq: A SeqRecord object corresponding to the target/reference sequence against which
            mutant sequences are compared.
        query_seqs: A list of SeqRecord objects corresponding to a collection of mutated sequences.
        kwargs: key word arguments to be forwarded to Bio.Align.Applications.MafftCommandline.

    Returns:
        a MultipleSeqAlignment object, where the first line corresponds to the target/reference
        sequence.
    """
    with tempfile.NamedTemporaryFile() as temp:
        SeqIO.write([target_seq] + query_seqs, temp.name, "fasta")
        mafft_cline = MafftCommandline(input=temp.name, **kwargs)
        stdout, stderr = mafft_cline()

    return AlignIO.read(io.StringIO(stdout), "fasta", alphabet=target_seq.seq.alphabet)


def get_mutations(
    target_seq_ali: SeqRecord, query_seq_ali: SeqRecord, offset: int = DEFAULT_OFFSET
) -> List[Mutation]:
    """
    Determine mutations between aligned target and query sequence.

    Args:
        target_seq_ali: A SeqRecord object corresponding to an aligned target/reference sequence.
        query_seq_ali: A SeqRecord object corresponding to an aligned mutant sequence.
        offset: Positional offset to use when describing mutations' position.

    Returns:
        A list of Mutation objects, describing each position in the alignment where the target and
        query sequences disagree.

    Raises:
        ValueError if the two sequences have different lengths.
    """
    if len(target_seq_ali) != len(query_seq_ali):
        raise ValueError("Sequences must be the same length")

    mutations = []

    i_target = 0
    i_query = 0

    # Identify mutations by iterating over aligned sequences.
    for i, (target, query) in enumerate(
        zip(str(target_seq_ali.seq), str(query_seq_ali.seq))
    ):
        has_mutation = target.lower() != query.lower()
        if has_mutation:
            mutations.append(
                Mutation(
                    target_seq_ali.id,
                    query_seq_ali.id,
                    i,
                    i_target,
                    i_query,
                    target,
                    query,
                    target + str(i_target) + query,
                    0,
                )
            )
        # Increment position counters.
        if target is GAP_SYM:
            i_query += 1
        elif query is GAP_SYM:
            i_target += 1
        else:
            i_target += 1
            i_query += 1

    if offset != 0:
        mutations = shift_mutations(mutations, offset)

    return mutations


def shift_mutations(mutations: List[Mutation], shift: int) -> List[Mutation]:
    """
    Shift the positions of a list of mutations.

    Arg:
        mutations: a list of Mutations.
        shift: the amount to shift the positions.

    Returns:
        a list of Mutations, with positions shifted.
    """
    return [
        Mutation(
            m.target_id,
            m.query_id,
            m.alignment_position + shift,
            m.target_position + shift,
            m.query_position + shift,
            m.target_seq,
            m.query_seq,
            m.target_seq + str(m.target_position + shift) + m.query_seq,
            m.offset + shift,
        )
        for m in mutations
    ]


def get_mutations_offset(mutations: List[Mutation]) -> int:
    """
    Get the offset value for a list of Mutations.

    Args:
        mutations: a list of Mutations.

    Returns:
         the offset of the Mutations.

    Raises:
        ValueError if the offsets aren't all the same.
    """
    offsets = set(m.offset for m in mutations)
    if len(offsets) != 1:
        raise ValueError("All mutations must have the same offset")
    return offsets.pop()


def reset_mutations_offset(mutations: List[Mutation]) -> List[Mutation]:
    """
    Reset  Mutations' position offset to zero.

    Args:
        mutations: a list of Mutations.

    Returns:
        a list of Mutations, with positions shifted back to zero.
    """
    return shift_mutations(mutations, -get_mutations_offset(mutations))


def msa_to_sequence_df(msa: MultipleSeqAlignment) -> pd.DataFrame:
    """
    Convert a multiple sequence alignment object to a pandas dataframe where each row corresponds to
    a sequence in the alignment and each column to a position in the alignment.

    Example::

        >>> print(msa)
        Alphabet() alignment with 3 rows and 5 columns
        actga s0
        acaaa s1
        acaat s2

        >>> df = msa_to_sequence_df(msa)
        >>> print(df)
            0  1  2  3  4
        s0  a  c  t  g  a
        s1  a  c  a  a  a
        s2  a  c  a  a  t

    Args:
        msa: a Bio.Align.MultipleSeqAlignment object.

    Returns:
        an n-sequences by sequence-length DataFrame, indexed on sequence id.
    """
    return pd.DataFrame(data=[list(r.seq) for r in msa], index=[r.id for r in msa])


def msa_to_mutations_df(msa: MultipleSeqAlignment) -> pd.DataFrame:
    """
    Convert a multiple sequence alignment object to a pandas dataframe where each row corresponds
    to one mutation of a query sequence relative to the target sequence.

    This method assumes that the first sequence in the alignment is the target/reference.

    Example::

        >>> print(msa)
        Alphabet() alignment with 3 rows and 5 columns
        actga s0
        acaaa s1
        acaat s2

        >>> df = msa_to_sequence_df(msa)
        >>> print(df)
          target_id query_id  alignment_position  target_seq query_seq name ...
        0        s0       s1                   2           t         a  t2a
        1        s0       s1                   3           g         a  g3a
        2        s0       s2                   2           t         a  t2a
        3        s0       s2                   3           g         a  g3a
        4        s0       s2                   4           a         t  a4t

    Args:
        msa: a Bio.Align.MultipleSeqAlignment object.

    Returns:
        an n-mutations by m-features DataFrame

    """
    ref = msa[0]
    aligned_seqs = msa[1:]
    mutations = [
        mutation for seq in aligned_seqs for mutation in get_mutations(ref, seq)
    ]

    return pd.DataFrame(mutations)


def msa_to_mutation_composition_df(msa: MultipleSeqAlignment) -> pd.DataFrame:
    """
    Convert a multiple sequence alignment object to a pandas dataframe where
    mutations are one-hot encoded

    This method assumes that the first sequence in the alignment is the target/reference.

    Example::

        >>> print(msa)
        Alphabet() alignment with 3 rows and 5 columns
        actga s0
        acaaa s1
        acaat s2

        >>> df = msa_to_mutation_copmosition_df(msa)
        >>> print(df)
            t2a  g3a  a4t
        s0    0    0    0
        s1    1    1    0
        s2    1    1    1

    Args:
        msa: a Bio.Align.MultipleSeqAlignment object.

    Returns:
        an n-sequences by m-mutations DataFrame
    """
    mutations_df = msa_to_mutations_df(msa)
    mutations_df["i"] = 1
    df = pd.pivot_table(
        mutations_df, index="query_id", columns="name", values="i", fill_value=0
    )
    df.index.name = "id"

    # Add rows for any sequence identical to the reference, including the reference.
    msa_index = [r.id for r in msa]
    missing_rows = list(set(msa_index) - set(df.index))
    missing_df = pd.DataFrame(index=missing_rows)
    df = df.append(missing_df, sort=False)
    df = df.fillna(0).astype(int)
    df = df.reindex(index=msa_index)

    # Sort columns by (position, from, to)
    cols = df.columns.values
    sorted_cols = sort_mutations(cols)
    df = df.loc[:, sorted_cols]

    return df


def sort_mutations(mutations: List[str]) -> List[str]:
    """
    Sort a list of mutation descriptors of the form "<from><position><to>" by (position, from, to).

    Args:
        mutations: a list of mutation descriptors.

    Returns:
        the sorted list of mutations.
    """
    sort_on = [
        (int(v[1]), v[0], v[2]) for v in [COMP_COL_SORT_RE.split(m) for m in mutations]
    ]
    sorted_mutations = list(zip(*sorted(zip(sort_on, mutations))))[1]
    return sorted_mutations


def msa_to_coverage_df(msa: MultipleSeqAlignment, mincount: int = 0) -> pd.DataFrame:
    """
    Convert a multiple sequence alignment to a pandas dataframe describing the number of occurrences
    of each type of mutation at each position.

    This method assumes that the first sequence in the alignment is the target/reference.

    Example::

        >>> print(msa)
        Alphabet() alignment with 3 rows and 5 columns
        actga s0
        acaaa s1
        acaat s2

        >>> df = msa_to_coverage_df(msa)
        >>> print(df)
           a  c  g  t
        0  0  0  0  0
        1  0  0  0  0
        2  2  0  0  0
        3  2  0  0  0
        4  0  0  0  1

    Args:
        msa: a Bio.Align.MultipleSeqAlignment object.
        mincount: include positions that contain at least this many mutations.

    Returns:
        a DataFrame
    """
    mutations_df = msa_to_mutations_df(msa)
    offset = mutations_df.offset[0]

    # Count mutations at each position.
    c = (
        mutations_df.groupby(["alignment_position", "query_seq"])["name"]
        .count()
        .reset_index()
    )
    c.columns = ["position", "mutation", "count"]
    c["mutation"] = c["mutation"].str.upper()

    # Pivot to position x mutation.
    coverage_df = c.pivot("position", "mutation", "count").fillna(0)

    # Determine valid letters
    # Sample letters from first "few" sequences.
    letters = "".join(set(c for seq in msa[:5000] for c in str(seq.seq)))

    # Add cols for any missing bases
    all_cols = sorted(letters.upper())
    missing_cols = set(all_cols) - set(coverage_df.columns)
    for c in missing_cols:
        coverage_df[c] = 0
    coverage_df = coverage_df[all_cols]
    coverage_df = coverage_df.astype(int)

    # Add rows for any missing positions
    all_pos = range(offset, len(msa[0]) + offset)
    missing_pos = set(all_pos) - set(coverage_df.index)
    missing_pos_df = pd.DataFrame(0, index=list(missing_pos), columns=all_cols)
    coverage_df = pd.concat([coverage_df, missing_pos_df], axis=0)
    coverage_df.sort_index(inplace=True)

    # Add reference letter to position index for easy reading.
    coverage_df.index = [f"{a.upper()} {i}" for i, a in enumerate(msa[0], offset)]
    coverage_df.index.name = "position"

    # Optionally drop positions with less than mincount observations.
    if mincount > 0:
        drop_rows = coverage_df.index[coverage_df.sum(axis=1) < mincount]
        coverage_df.drop(drop_rows, axis=0, inplace=True)

    return coverage_df


def msa_to_hamming_distance(msa: MultipleSeqAlignment) -> pd.DataFrame:
    """
    Calculate the hamming distance between all pairs of sequences in the msa.

    Args:
        msa: a Bio.Align.MultipleSeqAlignment object.

    Returns:
        a n-sequence by n-sequence DataFrame.
    """
    # Dist function requires numeric data :(
    seq_df = msa_to_sequence_df(msa).applymap(ord)

    # Compute hamming distance. Don't normalize by sequence length.
    hamm = dist.cdist(seq_df, seq_df, "hamming") * len(seq_df.iloc[0])
    hamm_df = pd.DataFrame(hamm, index=seq_df.index, columns=seq_df.index)
    hamm_df = hamm_df.astype(int)

    return hamm_df


def reverse_translate_mutations(
    mutations: List[Mutation], dna_sequences: List[SeqRecord], offset: int = 1
) -> List[Mutation]:
    """
    Convert amino acid mutation objects into a codon mutation objects.

    Note that reverse translation is based on actual dna sequence usage in provided
    sequences (as opposed to codon table look-ups).

    Args:
        mutations: A list of amino acid Mutation objects describing differences between aligned
            target and query amino acid sequences.
        dna_sequences: A list of SeqRecord objects whose ids match those referenced by the provided
            mutation objects
        offset: Positional offset to use when describing mutations' position.


    Returns:
        A list of Mutation objects, describing each input mutation in terms of codons.

    Raises:
        ValueError if mutations consist only consist of DNA.
        ValueError if dna_sequences do not consist of DNA.
    """
    # Ensure that mutations use 0-indexed.
    mutations = reset_mutations_offset(mutations)

    # TODO(flash): Add validation further up in the code to ensure sequence ids are unique.
    dna_sequences_map = {s.id: s for s in dna_sequences}

    # TODO(flash): These checks _should_ work for any real world example. Still, it would be better
    # to use Bio.Alphabet.
    # TODO(flash): Add sequence alphabet inference methods to avoid these heuristics.

    # Validate mutations are AAs.
    mutations_letters = "".join(
        _seq_to_str(l).upper() for m in mutations for l in (m.target_seq, m.query_seq)
    )
    if len(set(mutations_letters) - set(IUPACUnambiguousDNA.letters)) == 0:
        raise ValueError("mutations contain invalid letters.")

    # Validate dna_sequences are DNA.
    for seq in dna_sequences[:100]:
        if len(set(_seq_to_str(seq).upper()) - set(IUPACUnambiguousDNA.letters)) > 0:
            raise ValueError("dna sequences contain invalid letters.")

    codon_mutations = []
    for m in mutations:
        target_sequence = dna_sequences_map[m.target_id]
        target_start = m.target_position * 3
        target_end = target_start + 3 * len(m.target_seq)
        target_seq = _seq_to_str(target_sequence[target_start:target_end])
        query_sequence = dna_sequences_map[m.query_id]
        query_start = m.query_position * 3
        query_end = query_start + 3 * len(m.query_seq)
        query_seq = _seq_to_str(query_sequence[query_start:query_end])

        cm = Mutation(
            m.target_id,
            m.query_id,
            m.alignment_position * 3,
            target_start,
            query_start,
            target_seq,
            query_seq,
            target_seq + str(target_start) + query_seq,
            0,
        )

        codon_mutations.append(cm)

    if offset != 0:
        shift_mutations(codon_mutations, offset)

    return codon_mutations


def reverse_translate_mutations_df(
    mutations_df: pd.DataFrame, dna_sequences: List[SeqRecord]
) -> pd.DataFrame:
    """
    Reverse translate amino acid mutations dataframe generated by the msa_to_mutations_df()
    method to codon mutations.

    Note that reverse translation is based on actual dna sequence usage in the associated
    sequences, as opposed to codon table look-ups.

    Args:
        mutations_df: a dataframe of amino acid mutations, generated by the msa_to_mutations_df
            method.
        dna_sequences: a list of SeqRecord objects containing dna sequences referenced
            (by sequence id) by the mutations df.

    Returns:
        a dataframe describing mutations in terms of codons.
    """
    mutations = [Mutation(*row[1:]) for row in mutations_df.itertuples()]
    codon_mutations = reverse_translate_mutations(mutations, dna_sequences)

    return pd.DataFrame(codon_mutations)


def get_codon_usage_df(
    mutations_df, mutation_col: str, codon_col: Union[str, List]
) -> pd.DataFrame:
    """
    Calculate codon usage statistics for distinct mutation in the given mutations dataframe.

    The mutations dataframe is assumed to have been generated by the msa_to_mutation_copmosition_df
    method and augmented with codon-level mutation data determined by the
    reverse_translate_mutations_df method.

    The resulting usage dataframe includes:
        * <codon_col>_freq:  codon frequencies (as a dict) per mutation.
        * <codon_col>_mode: most frequent codon per mutation. In case of ties, multiple values will
            appear in a list.

    Args:
        mutations_df: a dataframe of amino acid mutations, generated by the msa_to_mutations_df
            method, and augmented with codon level information.
        mutation_col: the name of the column with unique mutation identifiers.
        codon_col: the name of one or more codon columns to summarize.

    Returns:
        a dataframe of codon usage statistics.
    """
    # Coerce inputs col to be list-like.
    if isinstance(codon_col, str):
        codon_col = [codon_col]

    # Iterate over columns of interest.
    usage_data = []
    for col_name in codon_col:
        # Get most common value per mutation.
        s1 = mutations_df.groupby(mutation_col)[col_name].agg(pd.Series.mode)
        s1.name = f"{col_name}_mode"
        usage_data.append(s1)

        # Get value counts per mutation.
        s2 = mutations_df.groupby(mutation_col)[col_name].agg(
            lambda x: x.value_counts().to_dict()
        )
        s2.name = f"{col_name}_freq"
        usage_data.append(s2)

    # Concat stats data.
    return pd.concat(usage_data, axis=1)


def get_covariate_mutations(
    composition_df: pd.DataFrame, min_covariates: int = 2
) -> Dict[str, List[str]]:
    """
    Identify mutations (columns) which co-vary within the given composition dataframe.

    Args:
        composition_df: a dataframe of one-hot encoded sequence mutation data, for instance, as
            generated by the msa_to_mutation_composition_df method.
        min_covariates: only include columns for which the number of covariates is at least this
            many.

    Returns:
        a dict whose keys are the first co-varying column and values are the complete list
        of covariates.
    """
    # Drop row duplicates to avoid indexing issues.
    df = composition_df.drop_duplicates().T
    df.index.name = "index"

    # Generate mapping using pandas magic...
    mapping = (
        df.reset_index()
        .groupby(df.columns.tolist())["index"]
        .agg(["first", list])
        .set_index("first")["list"]
        .to_dict()
    )

    # Filter by min count.
    mapping = {c: cvs for c, cvs in mapping.items() if len(cvs) >= min_covariates}

    return mapping


def _seq_to_str(s: Union[str, Seq, SeqRecord]) -> str:
    """
    Convert seq-like objects to strings.
    """
    if isinstance(s, SeqRecord):
        return str(s.seq)
    elif isinstance(s, Seq):
        return str(s)
    else:
        return str(s)
