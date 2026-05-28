import unittest

from rnai_designer.design import DesignConfig, build_construct_sequence, design_candidates, reverse_complement
from rnai_designer.fasta import FastaRecord


class DesignTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")

    def test_build_construct_sequence(self):
        self.assertEqual(build_construct_sequence("AACCGG", "TT"), "AACCGGTTCCGGTT")

    def test_design_candidates_keeps_shared_targets_and_rejects_offtargets(self):
        shared = "ACGT" * 80
        hap1 = FastaRecord("Hap1_LBD4", "Hap1_LBD4", shared)
        hap2 = FastaRecord("Hap2_LBD4", "Hap2_LBD4", "TTTTT" + shared + "AAAAA")
        paralog = FastaRecord("Other_LBD", "Other_LBD", "GATTACA" * 50)

        candidates = design_candidates(
            targets=[hap1, hap2],
            transcriptome=[hap1, hap2, paralog],
            config=DesignConfig(
                min_len=80,
                max_len=90,
                step=10,
                min_shared_fraction=1.0,
                min_efficient_sirnas=1,
                max_offtarget_mismatches=0,
                max_offtarget_transcripts=0,
                max_candidates=3,
            ),
        )

        self.assertTrue(candidates)
        self.assertEqual(candidates[0].target_stats["Hap1_LBD4"].shared_sirnas, candidates[0].total_sirnas)
        self.assertEqual(candidates[0].target_stats["Hap2_LBD4"].shared_sirnas, candidates[0].total_sirnas)
        self.assertEqual(candidates[0].offtarget_transcript_count, 0)

    def test_design_candidates_filters_offtargets(self):
        shared = "ACGT" * 80
        hap1 = FastaRecord("Hap1_LBD4", "Hap1_LBD4", shared)
        hap2 = FastaRecord("Hap2_LBD4", "Hap2_LBD4", shared)
        paralog = FastaRecord("Other_LBD", "Other_LBD", shared)

        candidates = design_candidates(
            targets=[hap1, hap2],
            transcriptome=[hap1, hap2, paralog],
            config=DesignConfig(
                min_len=80,
                max_len=80,
                step=10,
                min_shared_fraction=1.0,
                min_efficient_sirnas=1,
                max_offtarget_mismatches=0,
                max_offtarget_transcripts=0,
            ),
        )

        self.assertEqual(candidates, [])


if __name__ == "__main__":
    unittest.main()
