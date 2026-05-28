import unittest

from rnai_designer.bowtie import parse_bowtie_output


class BowtieTests(unittest.TestCase):
    def test_parse_bowtie_output(self):
        hits = parse_bowtie_output(
            [
                "q1\t+\ttranscriptA\t4\tACGT\tIIII\t0\t",
                "q2\t-\ttranscriptB\t10\tACGA\tIIII\t0\t3:A>T",
            ]
        )

        self.assertEqual(hits[0].query_id, "q1")
        self.assertEqual(hits[0].position, 5)
        self.assertEqual(hits[0].mismatches, 0)
        self.assertEqual(hits[1].transcript_id, "transcriptB")
        self.assertEqual(hits[1].mismatches, 1)


if __name__ == "__main__":
    unittest.main()
