"""test_quadratic_isogeny_primes.py

To run these tests, enter the following, possibly as sudo:

sage test_quadratic_isogeny_primes.py

    ====================================================================

    This file is part of Quadratic Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait

    Quadratic Isogeny Primes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    The author can be reached at: barinder.s.banwait@gmail.com

    ====================================================================

"""
#
# pylint: disable=no-self-use


import pytest
from sage.all import Integer, QuadraticField

from isogeny_primes import get_isogeny_primes
from sage_code.common_utils import CLASS_NUMBER_ONE_DISCS, EC_Q_ISOGENY_PRIMES

AUX_PRIME_COUNT = 10

# The first five examples are shown in the table in the article by
# Gonzaléz, Lario, and Quer

# The final case comes from Box's determination of quadratic points
# on X_0(73). From his table, we find that D = -31 should yield a
# 73. The other values in his table have either been checked above,
# or are class number one Ds.

# All these tests together less then 2 minutes on Sage v9.4
@pytest.mark.parametrize(
    "D, extra_isogeny, potenial_isogenies",
    [
        (-127, 73, {31, 41, 61, 79, 127}),
        (5 * 577, 103, {577}),
        (-31159, 137, {23, 29, 41, 61, 79, 109, 157, 313, 31159}),
        (61 * 229 * 145757, 191, {31, 173, 229, 241, 145757}),
        (11 * 17 * 9011 * 23629, 311, {71, 9011, 23629}),
        (-31, 73, {31, 41}),
    ],
)
def test_from_literature(D, extra_isogeny, potenial_isogenies):
    K = QuadraticField(D)
    upperbound = potenial_isogenies.union(EC_Q_ISOGENY_PRIMES).union({extra_isogeny})
    superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
    assert set(superset).issuperset(EC_Q_ISOGENY_PRIMES)
    assert extra_isogeny in superset
    assert upperbound.issuperset(superset)


# Check that the code actually runs for several Ds
R = 10


@pytest.mark.parametrize("D", range(-R, R + 1))
def test_interval(D):
    if Integer(D).is_squarefree() and D != 1:
        K = QuadraticField(D)
        if not K.discriminant() in CLASS_NUMBER_ONE_DISCS:
            superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
            assert set(superset).issuperset(EC_Q_ISOGENY_PRIMES)
