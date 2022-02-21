"""test_isogeny_primes.py

Test the top level function

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
    _, _ = get_isogeny_primes(K, **TEST_SETTINGS)
