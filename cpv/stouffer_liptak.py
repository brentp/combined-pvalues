
import numpy as np
from scipy.stats import norm, chisqprob
from numpy.linalg import cholesky as chol
from numpy.linalg.linalg import LinAlgError
import sys
qnorm = norm.ppf
pnorm = norm.cdf

def stouffer_liptak(pvals, sigma=None, correction=False):
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
    qvals = qnorm(1.0 - pvals, loc=0, scale=1).reshape(L, 1)
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
            result["OK"] = True
        except LinAlgError, e:
            print >>sys.stderr, e
            # cant do the correction non-invertible
            sigma *= 0.95
            np.fill_diagonal(sigma, 0.99)
            result["OK"] = False

        # http://en.wikipedia.org/wiki/Fisher's_method#Relation_to_Stouffer.27s_Z-score_method
    if correction:
        denom = np.sqrt(np.power(sigma, 2).sum())
        Cp = qvals.sum() / denom
    if not correction:
        Cp = qvals.sum() / np.sqrt(len(qvals))

    # get the right tail.
    pstar = norm.sf(Cp)
    if np.isnan(pstar):
        print >>sys.stderr, "BAD:", pvals, sigma
        pstar = np.median(pvals)
        result["OK"] = True
    result.update({"C": Cp, "p": pstar})
    return result

def z_score_combine(pvals, sigma):
    L = len(pvals)
    pvals = np.array(pvals, dtype=np.float64)
    pvals[pvals == 1] = 1.0 - 9e-16
    z = np.mean(qnorm(1.0 - pvals, loc=0, scale=1))
    lower = np.tril(sigma, k=-1) > 0
    sz = 1.0/L * np.sqrt(L + 2 * np.tril(sigma, k=-1).sum())
    res = {'p': norm.sf(z/sz), 'OK': True}
    return res

from scipy.misc import comb
from math import log, factorial
def zaykin_truncated_independent(pvals, cutoff=0.05):
    tau = cutoff
    pvals = np.asarray(pvals)
    L = len(pvals)
    # what is K?
    k = (pvals <= tau).sum()

    w = np.product(pvals[pvals <= tau]) # ? no qvalue, just a cutoff
    print w

    p = 0
    for k in range(1, L):
        lhs = comb(L, k) * (1 - tau)**(L - k)

        rhs = 0
        for s in range(k):
            rhs += ((k * log(tau) - log(w))**s) / factorial(s) * int(w <= tau**k) + tau**k * int(w > tau**k)


        p += lhs * w * rhs

    return p

def fisherp(pvals):
    """ combined fisher probability without correction """
    s = -2 * np.sum(np.log(pvals))
    return chisqprob(s, 2 * len(pvals))

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
