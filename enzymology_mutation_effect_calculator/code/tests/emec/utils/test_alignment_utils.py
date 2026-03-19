from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import unittest

from emec.utils import alignment_utils


class TestAlignmentUtils(unittest.TestCase):
    def test_align_multiple_sequences_mafft(self):
        ref_seq = SeqRecord(seq=Seq("actga"), id="0")
        mut_seq = SeqRecord(seq=Seq("acaaa"), id="1")

        msa = alignment_utils.align_multiple_sequences_mafft(ref_seq, [mut_seq])
        self.assertEqual(2, len(msa))
        self.assertEqual(5, len(msa[0]))
        self.assertEqual(str(ref_seq.seq), str(msa[0].seq))
        self.assertEqual(str(mut_seq.seq), str(msa[1].seq))

    def test_get_mutations(self):
        ref_seq = SeqRecord(seq=Seq("actga"), id="0")
        mut_seq = SeqRecord(seq=Seq("acaaa"), id="1")

        mutations = alignment_utils.get_mutations(ref_seq, mut_seq)

        self.assertEqual(2, len(mutations))
        self.assertEqual(3, mutations[0].alignment_position)
        self.assertEqual("t", mutations[0].target_seq)
        self.assertEqual("a", mutations[0].query_seq)
        self.assertEqual(4, mutations[1].alignment_position)
        self.assertEqual("g", mutations[1].target_seq)
        self.assertEqual("a", mutations[1].query_seq)

    def test_shift_mutations(self):
        ref_seq = SeqRecord(seq=Seq("actga"), id="0")
        mut_seq = SeqRecord(seq=Seq("acaaa"), id="1")

        mutations = alignment_utils.get_mutations(ref_seq, mut_seq)
        mutations = alignment_utils.shift_mutations(mutations, -1)

        self.assertEqual(2, len(mutations))
        self.assertEqual(2, mutations[0].alignment_position)
        self.assertEqual("t", mutations[0].target_seq)
        self.assertEqual("a", mutations[0].query_seq)
        self.assertEqual(3, mutations[1].alignment_position)
        self.assertEqual("g", mutations[1].target_seq)
        self.assertEqual("a", mutations[1].query_seq)

    def test_msa_to_sequence_df(self):
        ref_seq = SeqRecord(seq=Seq("actga"), id="0")
        mut_seq = SeqRecord(seq=Seq("acaaa"), id="1")
        msa = alignment_utils.align_multiple_sequences_mafft(ref_seq, [mut_seq])

        df = alignment_utils.msa_to_sequence_df(msa)

        self.assertEqual(str(ref_seq.seq), "".join(df.iloc[0].values))
        self.assertEqual(str(mut_seq.seq), "".join(df.iloc[1].values))

    def test_msa_to_mutations_df(self):
        ref_seq = SeqRecord(seq=Seq("actga"), id="0")
        mut_seq = SeqRecord(seq=Seq("acaaa"), id="1")
        msa = alignment_utils.align_multiple_sequences_mafft(ref_seq, [mut_seq])

        df = alignment_utils.msa_to_mutations_df(msa)

        self.assertTupleEqual((2, 9), df.shape)
        self.assertListEqual(["t", "g"], list(df["target_seq"]))
        self.assertListEqual(["a", "a"], list(df["query_seq"]))

    def test_msa_to_mutation_composition_df(self):
        ref_seq = SeqRecord(seq=Seq("actga"), id="0")
        mut_seq = SeqRecord(seq=Seq("acaaa"), id="1")
        msa = alignment_utils.align_multiple_sequences_mafft(ref_seq, [mut_seq])

        df = alignment_utils.msa_to_mutation_composition_df(msa)

        self.assertTupleEqual((2, 2), df.shape)
        self.assertListEqual([0, 0], list(df.iloc[0]))
        self.assertEqual([1, 1], list(df.iloc[1]))

    def test_msa_to_coverage_df(self):
        ref_seq = SeqRecord(seq=Seq("actga"), id="0")
        mut_seq = SeqRecord(seq=Seq("acaaa"), id="1")
        msa = alignment_utils.align_multiple_sequences_mafft(ref_seq, [mut_seq])

        df = alignment_utils.msa_to_coverage_df(msa)
        self.assertTupleEqual((5, 4), df.shape)
        self.assertEqual(1, df.loc["T 3", "A"])
        self.assertEqual(1, df.loc["G 4", "A"])

    def test_msa_hamming_distance(self):
        seqs = [
            SeqRecord(seq=Seq("actga"), id="0"),
            SeqRecord(seq=Seq("accga"), id="1"),
            SeqRecord(seq=Seq("aacta"), id="2"),
        ]

        msa = MultipleSeqAlignment(seqs)

        df = alignment_utils.msa_to_hamming_distance(msa)

        self.assertEqual(1, df.loc["0", "1"])
        self.assertEqual(3, df.loc["0", "2"])
        self.assertEqual(2, df.loc["1", "2"])
        for i in df.index:
            self.assertEqual(0, df.loc[i, i])

    def test_reverse_translate_mutations(self):
        dna_seqs = [
            SeqRecord(seq=Seq("atggctgat"), id="0"),
            SeqRecord(seq=Seq("atggttgat"), id="1"),
        ]
        aa_seqs = [dna.translate(id=True) for dna in dna_seqs]

        msa = MultipleSeqAlignment(aa_seqs)
        mutations = alignment_utils.get_mutations(msa[0], msa[1])

        codon_mutations = alignment_utils.reverse_translate_mutations(
            mutations, dna_seqs
        )
        cm = codon_mutations[0]

        self.assertEqual(1, len(codon_mutations))
        self.assertEqual("0", cm.target_id)
        self.assertEqual("1", cm.query_id)
        self.assertEqual(3, cm.alignment_position)
        self.assertEqual(3, cm.target_position)
        self.assertEqual(3, cm.query_position)
        self.assertEqual("gct", cm.target_seq)
        self.assertEqual("gtt", cm.query_seq)

    def test_reverse_translate_mutations_invalid_dna_sequences(self):
        dna_seqs = [
            SeqRecord(seq=Seq("atggctgat"), id="0"),
            SeqRecord(seq=Seq("atggttgat"), id="1"),
        ]
        aa_seqs = [dna.translate(id=True) for dna in dna_seqs]

        msa = MultipleSeqAlignment(aa_seqs)
        mutations = alignment_utils.get_mutations(msa[0], msa[1])

        with self.assertRaises(ValueError):
            alignment_utils.reverse_translate_mutations(mutations, aa_seqs)

    def test_reverse_translate_mutations_invalid_mutations(self):
        seqs = [
            SeqRecord(seq=Seq("actga"), id="0"),
            SeqRecord(seq=Seq("accga"), id="1"),
            SeqRecord(seq=Seq("aacta"), id="2"),
        ]

        msa = MultipleSeqAlignment(seqs)
        mutations = alignment_utils.get_mutations(msa[0], msa[1])

        with self.assertRaises(ValueError):
            alignment_utils.reverse_translate_mutations(mutations, seqs)

    def test_reverse_translate_mutations_df(self):
        dna_seqs = [
            SeqRecord(seq=Seq("atggctgat"), id="0"),
            SeqRecord(seq=Seq("atggttgat"), id="1"),
        ]
        aa_seqs = [dna.translate(id=True) for dna in dna_seqs]

        msa = MultipleSeqAlignment(aa_seqs)
        mutations_df = alignment_utils.msa_to_mutations_df(msa)

        df = alignment_utils.reverse_translate_mutations_df(mutations_df, dna_seqs)

        self.assertTupleEqual((1, 9), df.shape)
        self.assertListEqual(["gct"], list(df["target_seq"]))
        self.assertListEqual(["gtt"], list(df["query_seq"]))
        self.assertListEqual([3], list(df["alignment_position"]))

    def test_get_codon_usage_df(self):
        # Setup
        dna_seqs = [
            SeqRecord(seq=Seq("atg" "ggg" "ttt"), id="0"),
            # AA mutation: G1V
            SeqRecord(seq=Seq("atg" "gtc" "ttt"), id="1"),
            SeqRecord(seq=Seq("atg" "gtc" "ttt"), id="2"),
            SeqRecord(seq=Seq("atg" "gta" "ttt"), id="3"),
            # AA mutation: G1K
            SeqRecord(seq=Seq("atg" "aaa" "ttt"), id="4"),
        ]
        aa_seqs = [dna.translate(id=True) for dna in dna_seqs]
        msa = MultipleSeqAlignment(aa_seqs)
        mutations_df = alignment_utils.msa_to_mutations_df(msa)
        codon_df = alignment_utils.reverse_translate_mutations_df(
            mutations_df, dna_seqs
        )
        mutations_df["query_seq_codon"] = codon_df["query_seq"]

        # Get codon usage df.
        codon_col = "query_seq_codon"
        df = alignment_utils.get_codon_usage_df(mutations_df, "name", codon_col)

        # Check the dataframe in general.
        self.assertListEqual(["G2K", "G2V"], df.index.tolist())
        self.assertListEqual(
            [codon_col + "_mode", codon_col + "_freq"], df.columns.tolist()
        )

        # Check specific mutation.
        mutation = "G2V"
        self.assertEqual("gtc", df.loc[mutation, codon_col + "_mode"])
        self.assertDictEqual(
            {"gtc": 2, "gta": 1}, df.loc[mutation, codon_col + "_freq"]
        )

        # Check specific mutation.
        mutation = "G2K"
        self.assertEqual("aaa", df.loc[mutation, codon_col + "_mode"])
        self.assertDictEqual({"aaa": 1}, df.loc[mutation, codon_col + "_freq"])

    def test_get_codon_usage_df_with_tie(self):
        # Setup
        dna_seqs = [
            SeqRecord(seq=Seq("atg" "ggg" "ttt"), id="0"),
            # AA mutation: G1V
            SeqRecord(seq=Seq("atg" "gtc" "ttt"), id="1"),
            SeqRecord(seq=Seq("atg" "gta" "ttt"), id="2"),
            # AA mutation: G1K
            SeqRecord(seq=Seq("atg" "aaa" "ttt"), id="4"),
        ]
        aa_seqs = [dna.translate(id=True) for dna in dna_seqs]
        msa = MultipleSeqAlignment(aa_seqs)
        mutations_df = alignment_utils.msa_to_mutations_df(msa)
        codon_df = alignment_utils.reverse_translate_mutations_df(
            mutations_df, dna_seqs
        )
        mutations_df["query_seq_codon"] = codon_df["query_seq"]

        # Get codon usage df.
        codon_col = "query_seq_codon"
        df = alignment_utils.get_codon_usage_df(mutations_df, "name", codon_col)

        # Check specific mutation.
        mutation = "G2V"
        self.assertListEqual(
            ["gta", "gtc"], df.loc[mutation, codon_col + "_mode"].tolist()
        )
        self.assertDictEqual(
            {"gtc": 1, "gta": 1}, df.loc[mutation, codon_col + "_freq"]
        )

    def test_get_covariate_mutations(self):
        ref_seq = SeqRecord(seq=Seq("actgat"), id="0")
        mut_seqs = [
            SeqRecord(seq=Seq("acaaat"), id="1"),
            SeqRecord(seq=Seq("acaaac"), id="2"),
            SeqRecord(seq=Seq("acaaag"), id="3"),
        ]
        msa = alignment_utils.align_multiple_sequences_mafft(ref_seq, mut_seqs)
        df = alignment_utils.msa_to_mutation_composition_df(msa)

        cov_mapping = alignment_utils.get_covariate_mutations(df, min_covariates=1)
        expected = {"t6g": ["t6g"], "t6c": ["t6c"], "t3a": ["t3a", "g4a"]}
        self.assertDictEqual(expected, cov_mapping)

        cov_mapping = alignment_utils.get_covariate_mutations(df, min_covariates=2)
        expected = {"t3a": ["t3a", "g4a"]}
        self.assertDictEqual(expected, cov_mapping)

        cov_mapping = alignment_utils.get_covariate_mutations(df, min_covariates=3)
        expected = dict()
        self.assertDictEqual(expected, cov_mapping)
