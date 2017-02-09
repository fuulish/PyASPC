import sys
import numpy as np
from scipy.special import binom
from fractions import gcd

class ASPC(object):
    def __init__(self, data, chainlnth=2, damp=None, debug=False):
        """
        initialize ASPC Python object
        """

        self.debug = debug

        self.chainlnth = chainlnth

        if damp is not None:
            self.damp = damp
        else:
            self.damp = self.get_default_damp(self.chainlnth)

        self.totlnth = chainlnth + 2

        self.coeffs = self.get_coefficients(chainlnth=self.chainlnth)

        self.history = []
        for i in range(self.totlnth):
            self.history.append([])

    def get_coefficients(self, chainlnth=None):
        """
        """

        if chainlnth is None:
            chainlnth = self.chainlnth

        coeffs = self.generate_coefficients(chainlnth)

        return coeffs

    def predict(self, data_history, coeffs):
        """
        """

        prdat = np.zeros_like(data_history[0])

        return prdat

    def correct(self, f):
        """
        """

        #call f to get correction on top of prediction

        crdat = []

        return crdat

    def get_final_solution(self, prdat, crdat, damp):
        """
        """

        return damp * crdat + (1.-damp) * prdat

    def generate_coefficients(self, lnth):
        """
        """

        coeffs = []

        ordpo = lnth + 1

        totlnth = self.totlnth

        for i in range(totlnth):
            k = i + 1
            coeffs.append(self.get_Bj(ordpo, k))

        return coeffs

    def get_Bj(self, n, k):
        """
        """

        nmrtr = k * binom(2 * n + 2, n + 1 - k)
        dnmntr = binom(2 * n, n)
        sgn = np.power(-1, k+1)

        #FUDO| checking against Baranyai/Kiss values (up to k = 4: http://pubs.acs.org/doi/full/10.1021/ct5009069)
        if self.debug:

            gcd_val = gcd(nmrtr, dnmntr)
            print("COEFF is %i %i / %i" %(sgn, nmrtr / gcd_val, dnmntr / gcd_val))

            frst = sgn * nmrtr / dnmntr
            scnd = np.power(-1, k+1) * k * binom (2 * n + 2, n + 1 - k ) / binom (2 * n, n);

            if frst != scnd:
                print("EEEEERRROR.....")
                sys.exit(1)

        return sgn * nmrtr / dnmntr
        #return np.power(-1, k+1) * k * binom (2 * n + 2, n + 1 - k ) / binom (2 * n, n);

    def get_default_damp(self, lnth):
        """
        """

        return (lnth + 2.) / (2.*lnth + 3.);

    def set_chainlnth(self, chainlnth):
        """
        """

        self.chainlnth = chainlnth
        self.totlnth = self.chainlnth + 2

        self.coeffs = self.get_coefficients(chainlnth=self.chainlnth)

        #FUDO| need to update all other involved quantities
        #FUDO| update history
