from minimal_representation import get_minimal_representation
import unittest


class MinimalRepresentationTests(unittest.TestCase): 
    """
    We'll this out with a bunch of examples of expected behavior
    """
    def test_simple_snv(self): 
        pos, ref, alt = 1001, 'A', 'T'
        expected = 1001, 'A', 'T'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)

    def test_multiallelic_indels(self):
        pos, ref, alt = 1001, 'CTCC', 'CCCC'
        expected = 1002, 'T', 'C'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'CTCC', 'CCC'
        expected = 1001, 'CT', 'C'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'CTCC', 'CTC'
        expected = 1002, 'TC', 'T'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'CTAG', 'CTG'
        expected = 1002, 'TA', 'T'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'CTCC', 'CTACC'
        expected = 1002, 'T', 'TA'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'TCAGCAGCAG', 'TCAGCAG'
        expected = 1001, 'TCAG', 'T'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'CTT', 'CTTT'
        expected = 1001, 'C', 'CT'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'CTT', 'C'
        expected = 1001, 'CTT', 'C'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'CTT', 'CT'
        expected = 1001, 'CT', 'C'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'AAAATATATATAT', 'A'
        expected = 1001, 'AAAATATATATAT', 'A'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'AAAATATATATAT', 'AATAT'
        expected = 1001, 'AAAATATAT', 'A'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)
        pos, ref, alt = 1001, 'ACACACACAC', 'AACAC'
        expected = 1001, 'ACACAC', 'A'
        self.assertEqual(get_minimal_representation(pos, ref, alt), expected)

if __name__ == '__main__': 
    unittest.main()
