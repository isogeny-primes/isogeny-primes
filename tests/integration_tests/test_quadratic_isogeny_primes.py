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

import unittest

from sage.all import Integer, QuadraticField

from isogeny_primes import get_isogeny_primes
from sage_code.common_utils import CLASS_NUMBER_ONE_DISCS, EC_Q_ISOGENY_PRIMES

AUX_PRIME_COUNT = 10


class TestQuadraticIsogenyPrimes(unittest.TestCase):

    # The first five examples are shown in the table in the article by
    # Gonzaléz, Lario, and Quer

    # All these tests together take a few minutes on Sage v9.4
    def test_73(self):
        K = QuadraticField(-127)
        superset = get_isogeny_primes(K, AUX_PRIME_COUNT, loop_curves=True)
        self.assertTrue(set(superset).issuperset(EC_Q_ISOGENY_PRIMES))
        self.assertIn(73, superset)

    def test_103(self):
        K = QuadraticField(5 * 577)
        superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
        self.assertTrue(set(superset).issuperset(EC_Q_ISOGENY_PRIMES))
        self.assertIn(103, superset)

    def test_137(self):
        K = QuadraticField(-31159)
        superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
        self.assertTrue(set(superset).issuperset(EC_Q_ISOGENY_PRIMES))
        self.assertIn(137, superset)

    def test_191(self):
        K = QuadraticField(61 * 229 * 145757)
        superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
        self.assertTrue(set(superset).issuperset(EC_Q_ISOGENY_PRIMES))
        self.assertIn(191, superset)

    def test_311(self):
        K = QuadraticField(11 * 17 * 9011 * 23629)
        superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
        self.assertTrue(set(superset).issuperset(EC_Q_ISOGENY_PRIMES))
        self.assertIn(311, superset)

    # The next test comes from Box's determination of quadratic points
    # on X_0(73). From his table, we find that D = -31 should yield a
    # 73. The other values in his table have either been checked above,
    # or are class number one Ds.

    def test_73_Box(self):
        K = QuadraticField(-31)
        superset = get_isogeny_primes(K, AUX_PRIME_COUNT)
        self.assertTrue(set(superset).issuperset(EC_Q_ISOGENY_PRIMES))
        self.assertIn(73, superset)

    # Check that the code actually runs for several Ds
    def test_interval(self):

        R = 10

        for D in range(-R, R + 1):  # -86 takes very long, otherwise -100 to 100 works
            if Integer(D).is_squarefree() and D != 1:
                K = QuadraticField(D)
                if not K.discriminant() in CLASS_NUMBER_ONE_DISCS:
                    get_isogeny_primes(K, AUX_PRIME_COUNT)


if __name__ == "__main__":
    unittest.main(buffer=True)
