
import aspc
import numpy as np

def corrfun(prdat):
    return prdat * 2

for i in range(5):
    data = np.random.random(100)
    a = aspc.ASPC(data, chainlnth=i, debug=True, correction=corrfun)
    #c = a.get_coefficients()
    c = a.coeffs
    print c

for i in range(5):
    print i
    d = a.next()
    print d
