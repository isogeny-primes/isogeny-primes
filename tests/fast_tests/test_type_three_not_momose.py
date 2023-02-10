"""test_type_three_not_momose.py

Tests for TypeThreeNotMomosePrimes.

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
from sage_code.type_three_not_momose import type_three_not_momose
from sage_code.generic import get_strong_type_3_epsilons
from sage.all import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen

TEST_SETTINGS = {
    "bound": 1000,
}

x = polygen(QQ)
test_cases = [
    (x**2 + 1, 1, 1),
    (x**3 + 3 * x**2 + 1, 1, 0),
    (x**4 - 9 * x**2 - 10 * x + 50, 5, 1),
    (x**4 - x**3 - x**2 - 2 * x + 4, 1, 2),
]


@pytest.mark.parametrize("f, class_number, strong_L_count", test_cases)
def test_type_three_not_momose(f, class_number, strong_L_count):
    K = NumberField(f, "a")
    Kgal = K.galois_closure("b")
    embeddings = K.embeddings(Kgal)

    strong_type_3_epsilons = get_strong_type_3_epsilons(K, embeddings)

    assert len(strong_type_3_epsilons) == 2 * strong_L_count

    type_3_not_momose_primes, _ = type_three_not_momose(K, embeddings, strong_type_3_epsilons)

    if class_number == 1:
        assert type_3_not_momose_primes == []
