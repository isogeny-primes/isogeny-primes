import pytest
from sage.all import (
    QuadraticField,
    EllipticCurve_from_j,
    GF,
)

from sage_code.frobenius_polynomials import (
    semi_stable_frobenius_polynomial,
    isogeny_character_values_12,
)

K = QuadraticField(-127, "D")
D = K.gen(0)
j = 20 * (3 * (-26670989 - 15471309 * D) / 2 ** 26) ** 3
# this is one of the curves from the Gonzaléz, Lario, and Quer article
E = EllipticCurve_from_j(j)


def test_semi_stable_frobenius_polynomial():

    # test that we indeed have a 73 isogeny mod p
    for p in 2, 3, 5, 7, 11, 19:
        for pp, e in (p * K).factor():
            f = semi_stable_frobenius_polynomial(E, pp)
            assert not f.change_ring(GF(73)).is_irreducible()


@pytest.mark.parametrize(
    "p, values",
    [
        (2, {1, 8}),
        (3, {1}),
        (7, {72}),
        (11, {13 ** 12 % 73, 57 ** 12 % 73}),
    ],
)
def test_isogeny_character_values_12(p, values):
    for pp, _ in (p * K).factor():
        assert set(isogeny_character_values_12(E, 73, pp)) == values
