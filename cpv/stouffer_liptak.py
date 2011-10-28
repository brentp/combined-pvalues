
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
    """
    L = len(pvals)
    pvals = np.array(pvals)
    qvals = qnorm(1 - pvals, loc=0, scale=1).reshape(L, 1)
    # dont do the correction unless sigma is specified.
    result = {"OK": True}
    if not sigma is None:
        try:
            C = chol(sigma)
            Cm1 = np.matrix(C).I # C^-1
            # qstar
            qvals = Cm1 * qvals
        except LinAlgError:
            # cant do the correction non-invertible
            result["OK"] = False
    Cp = qvals.sum() / (len(qvals)**0.5)
    # get the right tail.
    pstar = 1 - pnorm(Cp)
    result.update({"C": Cp, "p": pstar})
    return result

def fisher(ps, aacf, df):
    # aacf is off-diagonal. 
    # copied from kechris.
    k0 = float(len(ps))
    v = float(df)
    psi = -2.0 * np.log(ps).sum()
    # TODO: numexpr
    cov = ( 3.263 * aacf + 0.710 * aacf**2 + 0.027 * aacf**3 + 0.727 * (1 / v)
          + 0.327 * (aacf / v) - 0.768 *(aacf**2) / v - 0.331 * (aacf**3) / v
          ).sum()
    ev = 2 * k0
    vv = 4 * k0 + 2 * cov
    f = 2 * ev**2 / vv
    cv = vv / (2 * ev)
    return chisqprob(psi / cv, f)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
