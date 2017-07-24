import sys
import numpy as np
from scipy.special import binom
from fractions import gcd

#FUDO| make more methods have less arguments, only those where external interference is expected

class ASPC(object):
    def __init__(self, data, chainlnth=2, damp=None, debug=False, correction=None, corrargs=()):
        """
        initialize ASPC Python object
        """

        self.debug = debug
        self.chainlength = chainlnth
        self.correction = correction
        self.corrargs = corrargs

        if damp is not None:
            self.damp = damp
        else:
            self.damp = ASPC.default_damp(self.chainlength)

        #FUDO| this is crap
        self.history = []
        for i in range(self._totlength):
            self.history.append(np.zeros_like(data))

        self.update_history(self.history, data)
        self.countme = 1

    @property
    def coeffs(self):
        """
        """
        return self._coeffs

    @coeffs.setter
    def coeffs(self, value):
        """
        """

        raise AttributeError("can't set coefficients, automagically handled by setting chainlength")

    def predict(self, history, coeffs):
        """
        """

        prdat = np.zeros_like(history[0])

        for c, h in zip(coeffs, history):
            prdat += c*h

        return prdat

    def correct(self, prdat, func, corrargs=()):
        """
        """

        #call f to get correction on top of prediction

        if func is None:
            sys.stdout('Cannot do correction without corresponding function\n')
            sys.exit(1)

        crdat = func(prdat, *corrargs)

        return crdat

    def get_final_solution(self, prdat, crdat, damp=None):
        """
        """

        if damp is None:
            damp = self.damp

        return damp * crdat + (1.-damp) * prdat

    def update_history(self, history, data):
        """
        """

        #FUDO| is that fine w.r.t. garbage collection later on?
        #FUDO| i.e., what will happen to history[-1] after call history = update_history(history, data)

        newhist = [data]

        for h in history[:-1]:
            newhist.append(h)

        self.history = newhist

    def next(self):
        """
        """

        #FUDO| should we use also a different total chain length? What would happen then?
        if self.countme <= self._totlength:
            self.coeffs = self.get_coefficients(self.countme - 2)
            print 'COEFFS', self.coeffs

        #FUDO| should I pass everyting as argument again?

        prdat = self.predict(self.history, self.coeffs)

        #the predicted data should be used in the correction function according to what is appropriate there

        crdat = self.correct(prdat, self.correction, self.corrargs)

        data = self.get_final_solution(prdat, crdat)

        self.update_history(self.history, data)

        self.countme += 1

        return data

    @staticmethod
    def generate_coefficients(length, totlength):
        """
        """

        coeffs = []

        ordpo = length + 1

        #FUDO| may get _totlength as an argument as well
        #FUDO| check whether lnth+2 > _totlength (consistency)

        for i in range(totlength):
            k = i + 1
            coeffs.append(ASPC.calculate_Bj(ordpo, k))

        return coeffs

    @staticmethod
    def calculate_Bj(n, k):
        """
        """

        nmrtr = k * binom(2 * n + 2, n + 1 - k)
        dnmntr = binom(2 * n, n)
        sgn = np.power(-1, k+1)

        #FUDO| checking against Baranyai/Kiss values (up to k = 4: http://pubs.acs.org/doi/full/10.1021/ct5009069)
        #if self.debug:

        #    gcd_val = gcd(nmrtr, dnmntr)
        #    print("COEFF is %i %i / %i" %(sgn, nmrtr / gcd_val, dnmntr / gcd_val))

        #    frst = sgn * nmrtr / dnmntr
        #    scnd = np.power(-1, k+1) * k * binom (2 * n + 2, n + 1 - k ) / binom (2 * n, n);

        #    if frst != scnd:
        #        print("EEEEERRROR.....")
        #        sys.exit(1)

        return sgn * nmrtr / dnmntr
        #return np.power(-1, k+1) * k * binom (2 * n + 2, n + 1 - k ) / binom (2 * n, n);

    @staticmethod
    def default_damp(lnth):
        """
        """

        return (lnth + 2.) / (2.*lnth + 3.);

    @property
    def chainlength(self):
        """
        I am the chain length of the predictor
        """
        return self._chainlength

    @chainlength.setter
    def chainlength(self, chainlength):
        """
        recalculates and sets the coefficients for a given chainlength
        """

        self._chainlength = chainlength
        self._totlength = self._chainlength + 2

        self._coeffs = ASPC.generate_coefficients(self._chainlength, self._totlength)

        #FUDO| need to update all other involved quantities
        #FUDO| update history

    @chainlength.deleter
    def chainlength(self):
        """
        """

        del self._chainlength
        del self._totlength
        del self._coeffs
