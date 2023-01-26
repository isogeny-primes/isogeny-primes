"""test_frobenius_polynomials.py

Test for the code in frobenius_polynomials.py

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2022 Barinder S. Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    The authors can be reached at: barinder.s.banwait@gmail.com and
    maarten@mderickx.nl.

    ====================================================================
"""

import pytest
from sage.all import EllipticCurve_from_j, GF, QQ, QuadraticField
from sage.rings.polynomial.polynomial_ring import polygen
from sage.schemes.elliptic_curves.constructor import EllipticCurve

from sage_code.frobenius_polynomials import (
    semi_stable_frobenius_polynomial,
    isogeny_character_values_12,
)

K = QuadraticField(-127, "D")
D = K.gen(0)
j = 20 * (3 * (-26670989 - 15471309 * D) / 2**26) ** 3
# this is one of the curves from the Gonzaléz, Lario, and Quer article
E = EllipticCurve_from_j(j)


def test_semi_stable_frobenius_polynomial():

    # test that we indeed have a 73 isogeny mod p
    for p in 2, 3, 5, 7, 11, 19:
        for pp, e in (p * K).factor():
            f = semi_stable_frobenius_polynomial(E, pp)
            assert not f.change_ring(GF(73)).is_irreducible()


def test_semi_stable_frobenius_polynomial_t():
    # an example where the frobenius polynomial depends on the purely ramified extension
    # we make
    x = polygen(QQ)
    K = QQ.extension(x - 1, "one")
    E = EllipticCurve(K, [49, 343])
    assert E.discriminant() == -(2**4) * 31 * 7**6
    assert E.j_invariant() == K(2**8 * 3**3) / 31
    f1 = semi_stable_frobenius_polynomial(E, K * 7, 1)
    f2 = semi_stable_frobenius_polynomial(E, K * 7, -1)(x=-x)
    assert f1 == f2


@pytest.mark.parametrize(
    "p, values",
    [
        (2, {1, 8}),
        (3, {1}),
        (7, {72}),
        (11, {13**12 % 73, 57**12 % 73}),
    ],
)
def test_isogeny_character_values_12(p, values):
    for pp, _ in (p * K).factor():
        assert set(isogeny_character_values_12(E, 73, pp)) == values
