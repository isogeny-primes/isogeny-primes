"""test_type_one_primes.py

Type 1 tests.

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

from sage_code.type_one_primes import (
    is_formal_immersion,
    cached_bad_formal_immersion_data,
)


@pytest.mark.parametrize(
    "d,p,expected_result",
    [
        [3, 71, 1],
        [4, 71, 1],
        [5, 71, 1],
        [6, 71, 1],
        [7, 71, 0],
    ],
)
def test_is_formal_immersion(d, p, expected_result):
    assert is_formal_immersion(d, p) == expected_result


@pytest.mark.parametrize("d", range(1, 13))
def test_cached_bad_formal_immersion_data(d):
    cached_bad_formal_immersion_data(d)
