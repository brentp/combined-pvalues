
import numpy as np
from scipy.stats import norm, chisqprob
from numpy.linalg import cholesky as chol
from numpy.linalg.linalg import LinAlgError
qnorm = norm.ppf
pnorm = norm.cdf

def stouffer_liptak(pvals, sigma=None):
    """
    The stouffer_liptak correction.
    >>> stouffer_liptak([0.1, 0.2, 0.8, 0.12, 0.011])
    {'p': 0.0168..., 'C': 2.1228..., 'OK': True}

    >>> stouffer_liptak([0.5, 0.5, 0.5, 0.5, 0.5])
    {'p': 0.5, 'C': 0.0, 'OK': True}

    >>> stouffer_liptak([0.5, 0.1, 0.5, 0.5, 0.5])
    {'p': 0.28..., 'C': 0.57..., 'OK': True}

    >>> stouffer_liptak([0.5, 0.1, 0.1, 0.1, 0.5])
    {'p': 0.042..., 'C': 1.719..., 'OK': True}

    >>> stouffer_liptak([0.5], np.matrix([[1]]))
    {'p': 0.5...}
    """
    L = len(pvals)
    pvals = np.asarray(pvals)
    qvals = qnorm(1 - pvals, loc=0, scale=1).reshape(L, 1)
    # dont do the correction unless sigma is specified.
    result = {"OK": True}
    if not sigma is None:
        try:
            C = chol(sigma)
            Cm1 = np.asmatrix(C).I # C^-1
            # qstar
            qvals = Cm1 * qvals
        except LinAlgError:
            # cant do the correction non-invertible
            result["OK"] = False
        # http://en.wikipedia.org/wiki/Fisher's_method#Relation_to_Stouffer.27s_Z-score_method
    """
        denom = np.sqrt(np.power(sigma, 2).sum())
        Cp = qvals.sum() / denom
    else:
    """
    Cp = qvals.sum() / np.sqrt(len(qvals))

    # get the right tail.
    pstar = 1 - pnorm(Cp)
    result.update({"C": Cp, "p": pstar})
    return result

def fisherp(pvals):
    """ combined fisher probability without correction """
    s = -2 * np.sum(np.log(pvals))
    return chisqprob(s, 2 * len(pvals))

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
