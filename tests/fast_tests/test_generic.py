"""test_generic.py

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
from sage.arith.misc import is_squarefree
from sage.rings.number_field.number_field import QuadraticField

from sage_code.generic import contains_imaginary_quadratic_field

square_free_D = [D for D in range(-100, 100) if is_squarefree(D) and D != 1]


@pytest.mark.parametrize("D", square_free_D)
def test_contains_imaginary_quadratic_field(D):
    K = QuadraticField(D)
    result = contains_imaginary_quadratic_field(K)
    assert result == (D < 0)
