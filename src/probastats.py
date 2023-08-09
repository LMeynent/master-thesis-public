import numpy as np
import scipy as sp
import scipy.stats as st
from scipy.integrate import quad

def gamma_pdf(mean=None, k=None, theta=None, alpha=None, beta=None, sd=None, var=None, cv=None):
    if alpha is not None:
        k = alpha
    
    if cv is not None:
        k = 1. / cv**2

    if beta is not None:
        theta = 1. / beta

    if sd is not None:
        var = sd**2

    if mean is not None:
        if var is not None:
            theta = var / mean
            k = mean / theta

        if k is not None:
            theta = mean / k

    if (k is None) or (theta is None):
        raise ValueError('Invalid combination of parameters')
    
    return lambda x: st.gamma.pdf(x, k, scale=theta)


def discretise_pdf(pdf, n=100):
    res = np.zeros(n)
    res[0] = quad(pdf, 0, 1.5)[0]
    for i in range(1, n):
        res[i] = quad(pdf, i+0.5, i+1.5)[0]
    return res / np.sum(res)