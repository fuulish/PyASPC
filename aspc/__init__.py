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
        self.chainlnth = chainlnth
        self.correction = correction
        self.corrargs = corrargs

        if damp is not None:
            self.damp = damp
        else:
            self.damp = self.get_default_damp(self.chainlnth)

        self.totlnth = chainlnth + 2

        self.coeffs = self.get_coefficients(chainlnth=self.chainlnth)

        #FUDO| this is crap
        self.history = []
        for i in range(self.totlnth):
            self.history.append(np.zeros_like(data))

        self.history = self.update_history(self.history, data)
        self.countme = 1

    def get_coefficients(self, chainlnth=None):
        """
        """

        if chainlnth is None:
            chainlnth = self.chainlnth

        coeffs = self.generate_coefficients(chainlnth)

        return coeffs

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

        return newhist

    def next(self):
        """
        """

        #FUDO| should we use also a different total chain length? What would happen then?
        if self.countme <= self.totlnth:
            self.coeffs = self.get_coefficients(self.countme - 2)
            print 'COEFFS', self.coeffs

        #FUDO| should I pass everyting as argument again?

        prdat = self.predict(self.history, self.coeffs)

        #the predicted data should be used in the correction function according to what is appropriate there

        crdat = self.correct(prdat, self.correction, self.corrargs)

        data = self.get_final_solution(prdat, crdat)

        self.history = self.update_history(self.history, data)

        self.countme += 1

        return data

    def generate_coefficients(self, lnth):
        """
        """

        coeffs = []

        ordpo = lnth + 1

        #FUDO| may get totlnth as an argument as well
        #FUDO| check whether lnth+2 > totlnth (consistency)

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
