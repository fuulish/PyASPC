from aspc import ASPC
import numpy as np

import unittest

a = ASPC(np.ones(1), chainlnth=2)

class TestASPCcoeffs(unittest.TestCase):

    def test_chain_two(self):
        self.assertEqual(a.coeffs[0], 2.8)
        self.assertEqual(a.coeffs[1], -2.8)
        self.assertEqual(a.coeffs[2], 1.2)
        self.assertEqual(a.coeffs[3], -0.2)

    #def test_chain_three(self):
    #def test_chain_four(self):

#class TestASPCpredict(unittest.TestCase):

if __name__ == '__main__':
    unittest.main()
