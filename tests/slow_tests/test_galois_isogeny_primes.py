"""test_galois_isogeny_primes.py

To run these tests, enter the following:

make integrationtests

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
from sage.all import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen

from isogeny_primes import get_isogeny_primes
from sage_code.common_utils import EC_Q_ISOGENY_PRIMES

TEST_SETTINGS = {
    "norm_bound": 50,
    "bound": 1000,
    "appendix_bound": 0,
    "ice_filter": False,
    "repeat_bound": 6,
}

"""
testcases generated getting various galois number fields from the lmfdb
"""

x = polygen(QQ)
test_cases = [
    (
        x ** 4 - x ** 3 + x ** 2 - x + 1,
        set(),
        {23, 29, 31, 41, 47, 53, 61, 73, 97, 103},
    ),
    (
        x ** 5 - x ** 4 - 4 * x ** 3 + 3 * x ** 2 + 3 * x - 1,
        set(),
        {23, 29, 31, 41, 47, 53, 59, 61, 71, 73, 89, 97, 109, 131},
    ),
    (
        x ** 6 + 13 * x ** 4 + 26 * x ** 2 + 13,
        set(),
        set(
            [23, 29, 31, 41, 47, 53, 59, 61, 73, 79, 83, 97, 103, 109, 113, 131]
            + [157, 181, 193, 233, 239, 281, 311, 313, 337, 359, 467, 677, 1747]
        ),
    ),
    (
        x ** 6 - 2 * x ** 5 + 6 * x ** 4 + 22 * x ** 3 + 41 * x ** 2 + 48 * x + 36,
        set(),

        set([23, 29, 31, 41, 47, 53, 59, 61, 71, 73, 79, 83, 97, 103, 109, 113] + [181, 191, 233, 251, 433, 499, 643]),

    ),
    (
        x ** 6 - 3 * x ** 5 - 2 * x ** 4 + 9 * x ** 3 - 5 * x + 1,
        set(),
        set(
            [23, 29, 31, 41, 47, 53, 59, 61, 71, 73, 79, 83, 97, 103, 107, 109, 113]
            + [127, 131, 137, 139, 199, 283, 439]
        ),
    ),
    (

        (x ** 7 - x ** 6 - 12 * x ** 5 + 7 * x ** 4) + (28 * x ** 3 - 14 * x ** 2 - 9 * x - 1),

        set(),
        set(
            [23, 29, 31, 41, 47, 53, 59, 61, 71, 73, 79, 83, 89, 97, 101, 103, 107]
            + [109, 113, 127, 131, 137, 157, 173, 211, 233, 311, 479, 4493]
            + [5849, 19121, 22721]
        ),
    ),
]


@pytest.mark.parametrize("f, extra_isogenies, potenial_isogenies", test_cases)
def test_galois(f, extra_isogenies, potenial_isogenies):
    K = NumberField(f, "a")

    superset, _ = get_isogeny_primes(K, **TEST_SETTINGS)
    assert set(EC_Q_ISOGENY_PRIMES).difference(superset) == set()
    assert extra_isogenies.difference(superset) == set(), "Known isogenies are filtered out!"
    upperbound = potenial_isogenies.union(EC_Q_ISOGENY_PRIMES).union(extra_isogenies)
    unlisted_potential_isogenies = superset.difference(upperbound)
    assert unlisted_potential_isogenies == set(), "We got worse at filtering"
    assert set(upperbound.difference(superset)) == set(), "We got better at filtering"
