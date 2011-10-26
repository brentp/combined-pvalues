import unittest
import os
import os.path as op
import sys
HERE = op.abspath(op.dirname(__file__))
sys.path.insert(0, op.abspath(op.join(HERE, "..", "..")))

BED = op.join(HERE, "data", "pvals.bed")

class CPVTest(unittest.TestCase):
    bed = BED

class TestACF(CPVTest):
    ranges = [1, 50, 100, 150]

    def test_acf(self):
        from cpv import acf
        res = acf.acf((self.bed,), self.ranges, 5)
        self.assert_(isinstance(res, list))
        self.assertEqual(len(res), len(self.ranges) - 1)

        last_cor = 100
        for i, row in enumerate(res):
            cor = row[1][0]
            self.assertEqual(row[0][0], self.ranges[i])
            self.assertEqual(row[0][1], self.ranges[i + 1])
            # only know this is true for the test data.
            self.assert_(cor < last_cor)
            cor = last_cor

        res = acf.acf((self.bed,), self.ranges, 5, partial=True)
        self.assertEqual(len(res), len(self.ranges) - 1)





if __name__ == "__main__":
    unittest.main()
