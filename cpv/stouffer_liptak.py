
import numpy as np
from scipy.stats import norm
from numpy.linalg import cholesky as chol
from numpy.linalg.linalg import LinAlgError
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
    result = {"OK": True}
    if not sigma is None:
        try:
            C = chol(sigma)
            Cm1 = np.matrix(C).I # C^-1
            # qstar
            qvals = Cm1 * np.asmatrix(qvals)
        except LinAlgError:
            # cant do the correction non-invertible
            result["OK"] = False
    Cp = qvals.sum() / (len(qvals)**0.5)
    pstar = pnorm(Cp)
    result.update({"C": Cp, "p": pstar})
    return result

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
