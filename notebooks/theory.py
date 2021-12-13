import numpy as np

def γ_zpv(β, r):
    frac = ((1+np.sqrt(6./r))/
            (1-np.sqrt(6./r)))
    K = -2.*(1.+2./r)**-1 * frac * (1.+frac)**-2
    return K*β**2