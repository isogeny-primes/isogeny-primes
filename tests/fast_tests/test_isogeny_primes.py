import pytest
from isogeny_primes import get_isogeny_primes
from sage.all import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

TEST_SETTINGS = {
    "norm_bound": 50,
    "bound": 1000,
    "ice_filter": True,
    "appendix_bound": 100,
    "auto_stop_strategy": True,
}

R = PolynomialRing(QQ, "x")
test_cases = [[1361, 0, 1]]


@pytest.mark.parametrize("coeffs", test_cases)
def test_get_isogeny_primes(coeffs):
    f = R(coeffs)
    K = QQ.extension(f, "a")
    _ = get_isogeny_primes(K, **TEST_SETTINGS)
