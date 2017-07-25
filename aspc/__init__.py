import sys
import numpy as np
from scipy.special import binom
from collections import deque

#FUDO| make more methods have less arguments, only those where external interference is expected

class ASPC(object):
    def __init__(self, data, chainlength=2, damp=None, debug=False, correction=None, corrargs=(), gradualstart=True):
        """
        initialize ASPC Python object
        """

        #FUX history and damp to properties as well?

        self._totlength = chainlength + 2
        self._chainlength = chainlength
        self._gradualstart = gradualstart

        self.history = deque([])
        self.update_history(data)

        self.update_chain()

        self.debug = debug

        self.correction = correction
        self.corrargs = corrargs
        self.damp = damp

        self.countme = 1

    @property
    def damp(self):
        """
        """
        return self._damp

    @damp.setter
    def damp(self, value):
        if value < 0 or value is None:
            self._damp = ASPC.default_damp(self.chainlength)
        else:
            self._damp = value

    @damp.deleter
    def damp(self):
        """
        """
        del self._damp

    @property
    def totlength(self):
        """
        """
        return self._totlength

    @totlength.setter
    def totlength(self, value):
        """
        """

        raise AttributeError("can't set total length of chain, please set chainlength attribute")

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

    def predict(self):
        """
        """

        prdat = np.zeros_like(self.history[0])

        for c, h in zip(self.coeffs, self.history):
            prdat += c*h

        return prdat

    def correct(self, prdat=None):
        """
        """

        #FUX| in principle this could happen, although it doesn't make much sense
        if prdat is None:
            raise RuntimeError('Cannot do correction without prediction, the point is to use prediction as initial guess for the correction')

        #FUX| for now it's a safety measure, define a correction or GTFO
        if self.correction is None:
            raise RuntimeError('Cannot do correction without corresponding function')

        #if self.correction is None:
        #    return prdat

        crdat = self.correction(prdat, *self.corrargs)

        return crdat

    def get_final_solution(self, prdat, crdat):
        """
        """

        return self.damp * crdat + (1.-self.damp) * prdat

    def update_history(self, data):
        """
        """

        if len(self.history) != 0:
            self.history.pop()

        self.history.appendleft(data)

    def next(self):
        """
        """

        #FUDO| should we use also a different total chain length? What would happen then?

        if self.gradualstart and self.countme <= self._totlength:
            self.chainlength = self.countme - 2
            #print 'COEFFS', self.coeffs

        #FUDO| should I pass everyting as argument again?

        prdat = self.predict()

        #the predicted data should be used in the correction function according to what is appropriate there

        crdat = self.correct(prdat)

        data = self.get_final_solution(prdat, crdat)

        self.update_history(data)

        self.countme += 1

        return data

    @staticmethod
    def generate_coefficients(length, totlength):
        """
        """

        coeffs = []

        ordpo = length + 1

        #FUDO| may get _totlength as an argument as well
        #FUDO| check whether length+2 > _totlength (consistency)

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

        return sgn * nmrtr / dnmntr

    @staticmethod
    def default_damp(length):
        """
        """

        return (length + 2.) / (2.*length + 3.);

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

        totlength = chainlength + 2

        if not self.gradualstart:
            if totlength < self._totlength:
                self._totlength = totlength

        if totlength > self._totlength:
            self._totlength = totlength

        self.update_chain()

    @chainlength.deleter
    def chainlength(self):
        """
        """

        del self._chainlength

    @property
    def gradualstart(self):
        return self._gradualstart

    @gradualstart.setter
    def gradualstart(self, value):

        self._gradualstart = value
        self.chainlength = self._chainlength

    def update_chain(self):

        self._coeffs = ASPC.generate_coefficients(self._chainlength, self._totlength)

        while self._totlength > len(self.history):
            self.history.append(np.zeros_like(self.history[0]))

        while self._totlength < len(self.history):
            self.history.pop()
