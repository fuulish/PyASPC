from aspc import ASPC
import numpy as np

import unittest


class TestASPCcoeffs(unittest.TestCase):
    def test_chain_zero(self):
        a = ASPC(np.ones(1), chainlength=0)
        self.assertEqual(a.coeffs[0], 2.0)
        self.assertEqual(a.coeffs[1], -1.0)

    def test_chain_one(self):
        a = ASPC(np.ones(1), chainlength=1)
        self.assertEqual(a.coeffs[0], 2.5)
        self.assertEqual(a.coeffs[1], -2.0)
        self.assertEqual(a.coeffs[2],  0.5)

    def test_chain_two(self):
        a = ASPC(np.ones(1), chainlength=2)
        self.assertEqual(a.coeffs[0], 2.8)
        self.assertEqual(a.coeffs[1], -2.8)
        self.assertEqual(a.coeffs[2], 1.2)
        self.assertEqual(a.coeffs[3], -0.2)

    def test_chain_three(self):
        a = ASPC(np.ones(1), chainlength=3)
        self.assertEqual(a.coeffs[0], 3.0)
        self.assertEqual(a.coeffs[1], -24./7.)
        self.assertEqual(a.coeffs[2], 27./14.)
        self.assertEqual(a.coeffs[3], -4./7.)
        self.assertEqual(a.coeffs[4], 1./14.)

    def test_chain_four(self):
        a = ASPC(np.ones(1), chainlength=4)
        self.assertEqual(a.coeffs[0], 22./7.)
        self.assertEqual(a.coeffs[1], -55./14.)
        self.assertEqual(a.coeffs[2], 55./21.)
        self.assertEqual(a.coeffs[3], -22./21.)
        self.assertEqual(a.coeffs[4], 5./21.)
        self.assertEqual(a.coeffs[5], -1./42.)

class TestASPCpredict(unittest.TestCase):
    def test_predict_correct_simple(self):
        def corrfun(prdat):
            return prdat

        data = np.ones(100)

        a = ASPC(data, chainlength=10, correction=corrfun)

        for i in range(1000):
            d = a.next()

        #self.assertEqual(d, data)
        np.testing.assert_allclose(d, data)

    def test_predict_correct(self):
        def corrfun(prdat):
            return prdat * 2

        data = np.random.random(100)

        a = ASPC(data, chainlength=2, correction=corrfun, gradualstart=False)

        for i in range(len(a.history)):
            a.history[i] = np.random.random(100)

        prd = np.zeros_like(data)
        coeffs = [2.8, -2.8, 1.2, -0.2]

        for c, h in zip(coeffs, a.history):
            prd += c * h

        cor = corrfun(prd)

        prd_aspc = a.predict()
        cor_aspc = a.correct(prd_aspc)

        #self.assertEqual(cor, nxt)
        np.testing.assert_allclose(cor, cor_aspc)

        nxt = a.damp * cor + (1.-a.damp) * prd
        nxt_aspc = a.next()

        np.testing.assert_allclose(nxt, nxt_aspc)

if __name__ == '__main__':
    unittest.main()
