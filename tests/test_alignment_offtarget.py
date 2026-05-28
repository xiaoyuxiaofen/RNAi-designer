import unittest

from rnai_designer.alignment import needleman_wunsch, project_reference_window
from rnai_designer.fasta import FastaRecord
from rnai_designer.offtarget import hamming_distance, scan_offtargets


class AlignmentAndOffTargetTests(unittest.TestCase):
    def test_alignment_projects_window_across_insertion(self):
        alignment = needleman_wunsch("hap1", "AAAACCCCGGGG", "hap2", "AAAATTTCCCCGGGG")
        projection = project_reference_window(alignment, 4, 8)

        self.assertEqual(projection.sequence, "CCCC")
        self.assertEqual(projection.target_id, "hap2")
        self.assertAlmostEqual(projection.aligned_identity, 1.0)

    def test_offtarget_scan_finds_one_mismatch_hit(self):
        query = "ACGTACGTACGTACGTACGTA"
        mutated = "ACGTACGTACATACGTACGTA"
        hits = scan_offtargets(
            [query],
            [FastaRecord("target", "target", query), FastaRecord("paralog", "paralog", mutated)],
            target_ids={"target"},
            max_mismatches=1,
        )

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].transcript_id, "paralog")
        self.assertEqual(hits[0].mismatches, 1)

    def test_hamming_distance_stops_after_threshold(self):
        self.assertEqual(hamming_distance("AAAA", "TTTT", stop_after=2), 3)


if __name__ == "__main__":
    unittest.main()
