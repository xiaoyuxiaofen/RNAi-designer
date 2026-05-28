import unittest

from rnai_designer.bowtie2 import parse_bowtie2_sam


class Bowtie2Tests(unittest.TestCase):
    def test_parse_bowtie2_sam(self):
        hits = parse_bowtie2_sam(
            [
                "@SQ\tSN:tx1\tLN:100",
                "q1\t0\ttx1\t5\t255\t21M\t*\t0\t0\tACGT\tIIII\tNM:i:0",
                "q2\t16\ttx2\t9\t255\t21M\t*\t0\t0\tACGA\tIIII\tNM:i:1",
            ]
        )

        self.assertEqual(hits[0].query_id, "q1")
        self.assertEqual(hits[0].position, 5)
        self.assertEqual(hits[0].strand, "+")
        self.assertEqual(hits[1].strand, "-")
        self.assertEqual(hits[1].mismatches, 1)


if __name__ == "__main__":
    unittest.main()
