from functools import lru_cache
import logging
from typing import Any, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from sqlalchemy.engine import Engine

from drawbridge.models.avro.dna_component import DnaComponent
from drawbridge.models.sample import Sample
from drawbridge.models.strain import Strain


def get_changes_table_reference_strain_id(strain_ids: List[int], engine: Engine) -> int:
    """
    Get the reference strain id (earliest ancestor) associated with a collection strains.

    This method assumes that strains of interest are present in the Genetic Changes Table.

    Args:
        strain_ids: a list of strain ids to query.

    Returns:
        the reference strain id.

    Raises:
        an Exception a unique reference strain id is not found.
    """
    q = f"""
    SELECT DISTINCT reference_strain_id as id
    FROM derived_data.genetic_changes gc
    WHERE gc.target_strain_id IN {sql_in_clause(strain_ids)}
    """
    df = pd.read_sql(q, engine)

    if len(df) != 1:
        raise Exception(
            f"Unable to locate unique strain id for input strains. Found: {df.id.tolist()}"
        )

    return df.iloc[0].id


def get_changes_table_logical_parent_strain_id(
    strain_ids: List[int], engine: Engine
) -> int:
    """
    Get the logical parent strain id associated with a collection strains.

    This method assumes that strians of interest are present in the Genetic Changes Table.

    Args:
        strain_ids: a list of strain ids to query.

    Returns:
        the logical parent strain id.

    Raises:
        an Exception if a unique logical parent strain id is not found.
    """
    q = f"""
    SELECT DISTINCT gc.logical_parent_id as id
    FROM derived_data.genetic_changes_dense gc
    WHERE gc.target_strain_id IN {sql_in_clause(strain_ids)}
    """
    df = pd.read_sql(q, engine)

    if len(df) != 1:
        raise Exception(
            f"Unable to locate unique strain id for input strains. Found: {df.id.tolist()}"
        )

    return df.iloc[0].id


def get_changes_table_payload_sequences(
    strain_ids: List[int],
    engine: Engine,
    feature_term: str = None,
    annotation_name: str = None,
    deduplicate_strategy: str = None,
) -> List[SeqRecord]:
    """
    Get payload sequences for a collection of strains based on "last edit" payload information in
    the Genetic Changes Table.

    If a feature term or feature term and annotation name are specified, check that the paylaod
    sequence annotation matches those qualifiers.

    Args:
        strain_ids: a list of strain ids to query.
        engine: an sqlalchemy engine connected to the specified database.
        feature_term: optionally filter the payload sequence by annotation feature_term.
        annotation_name: optionally filter the payload sequence by annotation name.
        deduplicate_strategy: if multiple payloads exist for a given strain, optionally drop duplicates
            using the given strategy:
            - None: Kepp all duplicates.
            - first: Keep the first occurrence.
            - last: Keep the last occurrence.
            - False: Drop all duplicates.

    Returns:
        a list of SeqRecord objects corresponding to the payload sequences. SeqRecord ids are those
        of the associated strains.
    """
    q = """
    SELECT DISTINCT
        target_strain_id as strain_id,
        target_strain_name as strain_name,
        payload_id
    FROM derived_data.genetic_changes_dense gc
    WHERE gc.target_strain_id IN {}
    """.format(
        sql_in_clause(strain_ids)
    )

    # Run the query.
    df = pd.read_sql(q, engine)

    # Get payload sequences safely.
    df["sequence"] = df.apply(
        lambda row: _safely(
            get_sequence, row.payload_id, feature_term, annotation_name
        ),
        axis=1,
    )
    df.dropna(subset=["sequence"], inplace=True)

    # Optionally, remove any duplicates.
    if deduplicate_strategy is not None:
        # Get ids with duplicates.
        counts_df = df.groupby("strain_id").size().reset_index(name="counts")
        dup_strain_ids = counts_df[counts_df["counts"] > 1]["strain_id"].to_list()

        if len(dup_strain_ids) > 0:
            logging.warning(
                f"Duplicate sequence data found for strains:\n {dup_strain_ids}\n"
                f"Removing duplicate data using strategy: keep='{deduplicate_strategy}'."
            )
            df.drop_duplicates(
                subset=["strain_id"], keep=deduplicate_strategy, inplace=True
            )
        else:
            logging.warning("No duplicate sequence data found.")

    # Convert the dataframe to seq records.
    seq_records = convert_dataframe_to_seq_records(
        df, "strain_id", "strain_name", "sequence"
    )

    return seq_records


@lru_cache(maxsize=64)
def get_sequence(
    zid: int,
    feature_term: str = None,
    annotation_name: str = None,
    chromosome_name=None,
) -> str:
    """
    Get the designated sequence for a given ZId (Sample, Strain, or DnaComponent).

    If a feature term or annotation name are specified, find the sub-sequence of
    the given Dna Component with matching annotation.

    Args:
        zid: The sequence identifier, either Sample, Strain, or DnaComponent.
        feature_term: optionally, return the sub-sequence with a feature_term annotation.
        annotation_name: optionally, return the sub-sequence with this feature_term annotation.
        chromosome_name: if zid is a Strain id and if that strain has multiple chromosomes,
          specify which one to use.

    Returns:
        A DNA sequence corresponding to the annotated payload region.

    Raises:
        a ValueError if a unique sequence matching the input args can not be identified.
    """

    # Get the strain or dna component id from sample id if that's what was  provided.
    # TODO(flash): This technique is brittle.  Use ZID.find(..).itemType
    if str(zid).startswith("9"):
        sample = Sample.find(zid)
        if sample.strain:
            zid = sample.strain.id
        elif sample.dnaComponent:
            zid = sample.dnaComponent.id
        else:
            raise ValueError(
                f"Sample must be associated with either a strain or dna component."
            )

    # Get the dna component associated with a strain id if that's what was provided.
    if str(zid).startswith("7"):
        strain = Strain.find(zid)
        dna_component = None

        # Handle various failure modes in which a unique chromosome id can not be determined.
        if not strain.genomeComponents:
            raise ValueError(
                f"Strain {zid} does not have any associated genome components."
            )
        elif chromosome_name and chromosome_name not in strain.genomeComponents:
            raise ValueError(
                f"Invalid chromosome name for strain {zid}. "
                f"Please specify one of {strain.genomeComponents.keys()}"
            )
        elif (
            not chromosome_name
            and len(strain.genomeComponents) > 1
            and not feature_term
        ):
            raise ValueError(
                f"Insufficient inputs to identify unique chromosome for strain {zid}."
                f"Please specify a specific chromosome or feature term."
            )

        # Handle various modes in which a preferred chromosome id specified or inferred.
        if chromosome_name and chromosome_name in strain.genomeComponents:
            # Get the dna component by chromosome name.
            dna_component_id = strain.genomeComponents[chromosome_name]
            dna_component = DnaComponent.find(dna_component_id)
        elif len(strain.genomeComponents) == 1 and not chromosome_name:
            # Get the one and only dna component.
            chromosome_name = list(strain.genomeComponents)[0]
            dna_component_id = strain.genomeComponents[chromosome_name]
            dna_component = DnaComponent.find(dna_component_id)
        elif len(strain.genomeComponents) > 1 and not chromosome_name and feature_term:
            # Search for a valid chromosome that contains...
            # Let's be lazy here and just return the first valid chromosome !?!
            for chromosome_name, dna_component_id in strain.genomeComponents.items():
                try:
                    get_sequence(dna_component_id, feature_term, annotation_name)
                    dna_component = DnaComponent.find(dna_component_id)
                    break
                except ValueError:
                    continue

        if not dna_component:
            raise ValueError(
                f"Insufficient inputs to identify unique chromosome for strain {zid}."
                f"Please specify a specific chromosome or feature term."
            )

    # If the dna component was provided directly, just look it up.
    elif str(zid).startswith("13"):
        dna_component = DnaComponent.find(zid)

    # Nothing to be done in this case.
    else:
        raise ValueError(f"Unable to determine sequence for ZId {zid}")

    # Apply annotation filters.
    if feature_term or (feature_term and annotation_name):
        dna_component = _find_payload_by_annotation_term_and_name(
            dna_component, feature_term, annotation_name
        )

    return dna_component.sequence


def _find_payload_by_annotation_term_and_name(
    dna_component: DnaComponent, feature_term: str, annotation_name: str = None
) -> DnaComponent:
    """
    Finds the payload sequence from a given dna_component by annotation feature term and
    name.

    The returned DnaComponent takes the name and description of the annotated region.

    If the payload sequence is on the reverse strand, return the reverse complement.

    Args:
        dna_component: The DnaComponent to search.
        feature_term: The feature_term to look for.
        annotation_name: Optional annotation name to look for.

    Returns:
        A DnaComponent corresponding to the annotated payload region.

    Raises:
        ValueError if a unique payload sequence can not be identified
    """
    annotations = [
        a
        for a in dna_component.list_annotations(feature_term=feature_term)
        if (not annotation_name or a.annotationName == annotation_name)
    ]
    if len(annotations) != 1:
        raise ValueError(
            "Unable to find a unique annotation {}:{}".format(
                feature_term, annotation_name
            )
            + "for DnaComponent {}.".format(dna_component.id)
        )

    payload_annotation = annotations[0]
    payload_component = dna_component[
        payload_annotation.startPosition : payload_annotation.endPosition
    ]
    payload_component.componentName = payload_annotation.annotationName
    payload_component.description = payload_annotation.description

    if payload_annotation.direction == DnaComponent.REVERSE:
        payload_component = DnaComponent.reverse_complement(payload_component)

    return payload_component


def sql_in_clause(ids: List[Any]) -> str:
    """
    Formats a list of values as "(val1,val2,...)" for use with a SQL IN clause.
    """
    return "(" + ",".join(str(i) for i in ids) + ")"


def convert_dataframe_to_seq_records(
    input_dataframe: pd.DataFrame, id_column: str, name_column: str, seq_column: str
) -> List[SeqRecord]:
    """
    Given a DataFrame with an id column, name column, and sequence column, returns a list of
    SeqRecords.

    Args:
        input_dataframe: The DataFrame containing id, name, and seq columns.
        id_column: The name of the column containing the id.
        name_column: The name of the column containing the name.
        seq_column: The name of the column containing the sequence.

    Returns:
        A list of SeqRecords contained in the DataFrame.
    """

    def parse_row_to_seq_record(row):
        id_value = str(row[id_column])
        name_value = row[name_column]
        seq_value = row[seq_column]
        seq = Seq(seq_value)
        return SeqRecord(seq, id=id_value, description=name_value)

    seq_records = input_dataframe.apply(parse_row_to_seq_record, axis=1).tolist()

    return seq_records


def _safely(f, *args, **kwargs):
    """
    Call a function safely by catching exceptions.

    Args:
        f: the function
        args: function args
        kwargs: function kwargs

    Returns:
        the function output or None if an except is raised.
    """
    try:
        return f(*args, **kwargs)
    except Exception:
        pass
