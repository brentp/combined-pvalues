from __future__ import print_function
import sys

import numpy as np
from scipy.stats import norm
import scipy.stats as ss
from numpy.linalg import cholesky as chol
from numpy.linalg.linalg import LinAlgError
qnorm = norm.ppf
pnorm = norm.cdf

chisqprob = ss.distributions.chi2.sf


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
    pvals = np.array(pvals, dtype=np.float64)
    pvals[pvals == 1] = 1.0 - 9e-16
    qvals = norm.isf(pvals, loc=0, scale=1).reshape(L, 1)
    if any(np.isinf(qvals)):
        raise Exception("bad values: %s" % pvals[list(np.isinf(qvals))])

    # dont do the correction unless sigma is specified.
    result = {"OK": True}
    if not sigma is None:
        try:
            C = chol(sigma)
            Cm1 = np.asmatrix(C).I # C^-1
            # qstar
            qvals = Cm1 * qvals
        except LinAlgError as e:
            result["OK"] = False
            result = z_score_combine(pvals, sigma)
            return result

    Cp = qvals.sum() / np.sqrt(len(qvals))
    # get the right tail.
    pstar = norm.sf(Cp)
    if np.isnan(pstar):
        print("BAD:", pvals, sigma, file=sys.stderr)
        pstar = np.median(pvals)
        result["OK"] = True
    result.update({"C": Cp, "p": pstar})
    return result

def z_score_combine(pvals, sigma):
    L = len(pvals)
    pvals = np.array(pvals, dtype=np.float64)
    pvals[pvals == 1] = 1.0 - 9e-16
    z = np.mean(norm.isf(pvals, loc=0, scale=1))
    sz = 1.0 /L * np.sqrt(L + 2 * np.tril(sigma, k=-1).sum())
    res = {'p': norm.sf(z/sz), 'OK': True}
    return res

def fisherp(pvals):
    """ combined fisher probability without correction """
    s = -2 * np.sum(np.log(pvals))
    return chisqprob(s, 2 * len(pvals))

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
