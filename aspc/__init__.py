import sys
import numpy as np
from scipy.special import binom
from collections import deque

class ASPC(object):
    def __init__(self, data, chainlength=2, damp=None, debug=False, correction=None, corrargs=(), gradualstart=True):
        """
        initialize ASPC Python object

        args:
            data            - initial guess (numpy.ndarray)
            chainlength     - length of predictor chain (int)
            damp            - damping factor to join predictor and corrector (float)
            debug           - turn on/off debug output (bool)
            correction      - function that performs correction step (function)
            corrargs        - additional arguments for correction function (list)
            gradualstart    - increase predictor chain according to (initially unfilled) history

        notes:
            The correction function takes as first argument the prediction in form of a numpy.ndarray. All other arguments
            to the correction function have to be passed through the corrargs list.

        """

        self._totlength = chainlength + 2
        self._chainlength = chainlength
        self._gradualstart = gradualstart

        self.history = deque([])
        self.update_history(np.array(data))

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
        """
        set damping factor
        arg:
            value - damping factor, if < 0 or None, then default parameter generated
        """

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
        predict new initial guess for corrector based on history of previous steps
        """

        prdat = np.zeros_like(self.history[0])

        for c, h in zip(self.coeffs, self.history):
            prdat += c*h

        return prdat

    def correct(self, prdat=None):
        """
        corrector step, one step in the self-consistent solution of the problem
        args:
            prdat - initial guess from predictor
        """

        #FUX| in principle this could happen, although it doesn't make much sense
        if prdat is None:
            raise RuntimeError('Cannot do correction without prediction, the point is to use prediction as initial guess for the correction')

        #FUX| for now it's a safety measure, define a correction or GTFO
        if self.correction is None:
            raise RuntimeError('Cannot do correction without corresponding function')

        crdat = self.correction(prdat, *self.corrargs)

        return crdat

    def get_final_solution(self, prdat, crdat):
        """
        combine predictor and corrector to get final result
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
        perform one complete cycle in predictor/corrector scheme
        """

        if self.gradualstart and self.countme <= self._totlength:
            self.chainlength = self.countme - 2

        prdat = self.predict()

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

        for i in range(totlength):
            k = i + 1
            coeffs.append(ASPC.calculate_Bj(ordpo, k))

        return coeffs

    @staticmethod
    def calculate_Bj(n, k):
        """
        calculate actual numeric values of coefficients used in predictor chain
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
        args:
            chainlength - (int)
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
        """
        """

        return self._gradualstart

    @gradualstart.setter
    def gradualstart(self, value):
        """
        """

        self._gradualstart = value
        self.chainlength = self._chainlength

    def update_chain(self):
        """
        remove old data, add placeholder data; all to fit to predictor chain length
        """

        self._coeffs = ASPC.generate_coefficients(self._chainlength, self._totlength)

        while self._totlength > len(self.history):
            self.history.append(np.zeros_like(self.history[0]))

        while self._totlength < len(self.history):
            self.history.pop()
