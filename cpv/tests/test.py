import unittest
import os
import os.path as op
import sys
import numpy as np
HERE = op.abspath(op.dirname(__file__))
sys.path.insert(0, op.abspath(op.join(HERE, "..", "..")))
from cpv import stouffer_liptak as sl

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

class TestStoufferLiptak(CPVTest):
    pvals = np.array([0.2, 0.05, 0.2])

    def testSLOK(self):
        res = sl.stouffer_liptak(self.pvals)
        self.assert_(all(k in res for k in "OK p C".split()))

    def testSigma(self):
        unadjusted = sl.stouffer_liptak(self.pvals)
        sigma = np.eye(3)
        sigma[0,1] = sigma[1,0] = 0.9
        sigma[2,1] = sigma[1,2] = 0.7
        sigma[2,2] = sigma[1,2] = 0.5
        sigma[2,0] = sigma[0,2] = 0.6
        adjusted = sl.stouffer_liptak(self.pvals, sigma)
        self.assert_(adjusted['p'] > unadjusted['p'])


if __name__ == "__main__":
    unittest.main()
