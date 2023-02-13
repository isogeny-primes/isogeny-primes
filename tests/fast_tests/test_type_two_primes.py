"""test_type_two_primes.py

Tests for type two primes.

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
from sage_code.type_two_primes import type_2_primes, get_type_2_not_momose
from sage.all import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen

TEST_SETTINGS = {
    "bound": 1000,
}

x = polygen(QQ)
test_cases = [
    (x**2 + 1, 1),
    (x**3 + 3 * x**2 + 1, 1),
    (x**4 - 9 * x**2 - 10 * x + 50, 5),
]


@pytest.mark.parametrize("f, class_number", test_cases)
def test_type_2_primes(f, class_number):
    K = NumberField(f, "a")
    Kgal = K.galois_closure("b")
    embeddings = K.embeddings(Kgal)

    type_two_not_momose = get_type_2_not_momose(K, embeddings)
    _ = type_2_primes(K, embeddings, **TEST_SETTINGS)

    if class_number == 1:
        assert type_two_not_momose == set()
