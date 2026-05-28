import unittest

from rnai_designer.accessibility import heuristic_accessibility, parse_lunp


class AccessibilityTests(unittest.TestCase):
    def test_parse_lunp_selects_requested_unpaired_length(self):
        text = """#i$ l=1 l=2 l=3
1 0.10 0.20 0.30
2 0.11 0.21 0.31
"""
        self.assertEqual(parse_lunp(text, 3), {1: 0.30, 2: 0.31})

    def test_heuristic_accessibility_bounds_values(self):
        values = heuristic_accessibility("ATGCATGC", 4)
        self.assertTrue(values)
        self.assertTrue(all(0.0 <= value <= 1.0 for value in values.values()))


if __name__ == "__main__":
    unittest.main()
