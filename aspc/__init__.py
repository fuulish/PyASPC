import numpy as np
from scipy.special import binom

class ASPC(object):
    def __init__(self, data, chainlnth=2, damp=None):
        """
        initialize ASPC Python object
        """

        self.chainlnth = chainlnth

        if damp is not None:
            self.damp = damp
        else:
            self.damp = self.get_default_damp(self.chainlnth)

        self.totlnth = chainlnth + 2

        self.coeffs = self.get_coefficients()

        self.history = []
        for i in range(self.totlnth):
            self.history.append([])

    def get_coefficients(self):
        """
        """

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

        #call f to get prediction

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

        return np.power(-1, k+1) * k * binom (2 * n + 2, n + 1 - k ) / binom (2 * n, n);

    def get_default_damp(self, lnth):
        """
        """

        return (lnth + 2.) / (2.*lnth + 3.);

    def set_chainlnth(self, chainlnth):
        """
        """

        self.chainlnth = chainlnth
        self.totlnth = self.chainlnth + 2

        #FUDO| need to update all other involved quantities
        #FUDO| update history
