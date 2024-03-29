"""test_cubic_isogeny_primes.py

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
    "appendix_bound": 80,
    "ice_filter": True,
}

"""
testcases generated by:
[
    (
        pari(hilbert_class_polynomial(D * f ** 2)).polredabs(),
        set(D.prime_divisors()),
        set(),
    )
    for D, f in cm_orders(3)
    if f == 1
]
"""

x = polygen(QQ)
test_cases = [
    (x**3 - x**2 + 1, {23}, set()),
    (x**3 + x - 1, {31}, set()),
    (x**3 + 2 * x - 1, {59}, set()),
    (x**3 - x**2 + x - 2, {83}, set()),
    (x**3 - x**2 + 3 * x - 2, {107}, set()),
    (x**3 - x**2 + x + 2, {139}, set()),
    (x**3 - 2 * x - 3, {211}, set()),
    (x**3 + 4 * x - 1, {283}, set()),
    (x**3 - x**2 + 3 * x + 2, {307}, set()),
    (x**3 - x**2 + 3 * x - 4, {331}, set()),
    (x**3 - x**2 + x - 4, {379}, set()),
    (x**3 + 4 * x - 3, {499}, set()),
    (x**3 - x**2 - 3 * x - 4, {547}, set()),
    (x**3 - 2 * x - 5, {643}, set()),
    (x**3 - x**2 + 2 * x - 12, {883}, set()),
    (x**3 - x**2 - 7 * x + 12, {907}, set()),
]

bad_cubic_formal_immersion_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 43, 73}


@pytest.mark.parametrize("f, extra_isogenies, potenial_isogenies", test_cases)
def test_cm_type_2(f, extra_isogenies, potenial_isogenies):
    K = NumberField(f, "a")

    superset, _ = get_isogeny_primes(K, **TEST_SETTINGS)
    assert set(EC_Q_ISOGENY_PRIMES).difference(superset) == set()
    # for p in extra_isogenies:
    #     assert satisfies_condition_CC(K, p), "Not of type 2"
    pnip = sorted(superset.difference(set(EC_Q_ISOGENY_PRIMES)))
    print("f = {}  disc = {} pnip = {}".format(f, K.discriminant(), pnip))
    assert extra_isogenies.difference(superset) == set(), "Known CM isogenies are filtered out!"
    upperbound = potenial_isogenies.union(EC_Q_ISOGENY_PRIMES).union(extra_isogenies)
    upperbound = upperbound.union(bad_cubic_formal_immersion_primes)
    unlisted_potential_isogenies = superset.difference(upperbound)
    assert len(unlisted_potential_isogenies) <= 2, "We got worse at filtering"
    if unlisted_potential_isogenies:
        assert max(unlisted_potential_isogenies) <= 109, "We got worse at filtering"
    # assert unlisted_potential_isogenies == set(), "We got worse at filtering"
    # assert set(upperbound.difference(superset)) == set(), "We got better at filtering"
