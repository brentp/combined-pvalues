
import numpy as np
from scipy.stats import norm
from numpy.linalg import cholesky as chol
qnorm = norm.ppf
pnorm = norm.cdf

def stouffer_liptak(pvals, sigma=None):
    """
    The stouffer_liptak correction.
    >>> stouffer_liptak([0.1, 0.2, 0.8,0.12, 0.11])
    {'p': 0.0497..., 'C': -1.647...}
    """
    L = len(pvals)
    qvals = qnorm(pvals, loc=0, scale=1).reshape(L, 1)
    # dont do the correction unless sigma is specified.
    if not sigma is None:
        C = chol(sigma)
        Cm1 = np.matrix(C).I # C^-1
        # qstar
        qvals = Cm1 * np.asmatrix(qvals)
    Cp = qvals.sum() / (len(qvals)**0.5)
    pstar = pnorm(Cp)
    return {"C": Cp, "p": pstar}

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
