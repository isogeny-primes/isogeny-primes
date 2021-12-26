from isogeny_primes import get_isogeny_primes
from .common_utils import EC_Q_ISOGENY_PRIMES
from sage.all import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

R = PolynomialRing(QQ, "x")
LMFDB_DEFAULTS = {
    "norm_bound": 50,
    "bound": 1000,
    "loop_curves": False,
    "use_PIL": False,
    "heavy_filter": True,
    "appendix_bound": 1000,
    "stop_strategy": "auto",
}


def isogeny_primes(coeffs, **kwargs):
    del kwargs["label"]
    f = R(coeffs)
    K = QQ.extension(f, "a")
    kwargs = {**LMFDB_DEFAULTS, **kwargs}
    primes = get_isogeny_primes(K, **kwargs)
    return sorted(EC_Q_ISOGENY_PRIMES) + sorted(primes.difference(EC_Q_ISOGENY_PRIMES))
