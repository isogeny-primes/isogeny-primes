"""test_weeding.py

Tests for weeding functions.

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

from sage.all import (
    QuadraticField,
    companion_matrix,
    prime_range,
)  # pylint: disable=no-name-in-module
from sage_code.weeding import (
    get_dirichlet_character,
    works_method_of_appendix,
    is_rank_of_twist_zero,
)


def test_get_dirichlet_character():
    K = QuadraticField(-31)
    chi = get_dirichlet_character(K)
    assert (3 * K).is_prime()
    assert chi(3) == -1
    assert not (5 * K).is_prime()
    assert chi(5) == 1
    assert (73 * K).is_prime()
    assert chi(73) == -1


@pytest.mark.parametrize(
    "D, p, works",
    [
        (5, 79, True),
        (-31, 73, False),
    ],
)
def test_works_method_of_appendix(D, p, works):
    K = QuadraticField(D)
    result = works_method_of_appendix(p, K)
    assert result == works


def test_is_rank_of_twist_zero():
    p = 73
    K = QuadraticField(-31)
    chi = get_dirichlet_character(K)
    assert not is_rank_of_twist_zero(p, chi)
