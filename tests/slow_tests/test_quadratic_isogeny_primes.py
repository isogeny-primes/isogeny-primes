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

# total running time of all tests in this file is about 5 minutes

# The first case comes from Box's determination of quadratic points
# on X_0(73). From his table, we find that D = -31 should yield a
# 73 isogeny. The other values in his table have either been checked
# Gonzaléz, Lario, and Quer or are class number one Ds.

# The last five examples are shown in the table in the article by
# Gonzaléz, Lario, and Quer

# The running time and the memory consumption of the method of the appendix is
# linear in the size of the conductor of the number field. Hence we
# set "appendix_bound" to 0 in order to disable the method of the appendix for
#   D = 61 * 229 * 145757 = 2036079533
#   D = 11 * 17 * 9011 * 23629 = 39816211853
# since in those cases just simply looping over all integers up to that
# bound is already taking ages, and even worse. Just list(range(D)) will exaust
# all memory on most machines.
@pytest.mark.parametrize(
    "D, extra_isogeny, appendix_bound, potenial_isogenies",
    [
        (-31, 73, 1000, {31, 41}),
        (-127, 73, 1000, {31, 41, 61, 79}),
        (5 * 577, 103, 1000, {577}),
        (-31159, 137, 1000, {23, 29, 41, 61, 79, 109, 31159}),
        (61 * 229 * 145757, 191, 0, {31, 173, 229, 241, 145757}),
        (11 * 17 * 9011 * 23629, 311, 0, {71, 9011, 23629}),
    ],
)
def test_from_literature(D, extra_isogeny, appendix_bound, potenial_isogenies):
    K = QuadraticField(D)
    upperbound = potenial_isogenies.union(EC_Q_ISOGENY_PRIMES).union({extra_isogeny})
    superset = get_isogeny_primes(K, AUX_PRIME_COUNT, appendix_bound=appendix_bound)
    assert set(superset).issuperset(EC_Q_ISOGENY_PRIMES)
    assert extra_isogeny in superset
    assert upperbound.issuperset(superset)
    assert set(upperbound.difference(superset)) == set()


# Check that the code actually runs for several Ds
R = 10


@pytest.mark.parametrize("D", range(-R, R + 1))
def test_interval(D):
    if Integer(D).is_squarefree() and D != 1:
        K = QuadraticField(D)
        if not K.discriminant() in CLASS_NUMBER_ONE_DISCS:
            superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
            assert set(superset).issuperset(EC_Q_ISOGENY_PRIMES)
