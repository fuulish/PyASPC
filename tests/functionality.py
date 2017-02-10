
import aspc
import numpy as np

debug = False

def corrfun(prdat):
    return prdat * 2

for i in range(5):
    data = np.random.random(100)
    a = aspc.ASPC(data, chainlnth=i, debug=debug, correction=corrfun)
    #c = a.get_coefficients()
    c = a.coeffs

data = np.linspace(0, 10, 100)
a = aspc.ASPC(data, chainlnth=2, debug=debug, correction=corrfun)
print a.history
print a.coeffs

for i in range(6):
    d = a.next()
    print d
