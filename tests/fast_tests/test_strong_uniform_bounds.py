"""test_strong_uniform_bounds.py

Test strong uniform bounds.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2023 Barinder S. Banwait and Maarten Derickx

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
from sage_code.strong_uniform_bounds import splitting_types


test_cases = [
    (
        5,
        [
            {"r": 1, "es": (1,), "fs": (5,)},
            {"r": 1, "es": (5,), "fs": (1,)},
            {"r": 2, "es": (1, 1), "fs": (4, 1)},
            {"r": 2, "es": (2, 1), "fs": (2, 1)},
            {"r": 2, "es": (4, 1), "fs": (1, 1)},
            {"r": 2, "es": (1, 1), "fs": (3, 2)},
            {"r": 2, "es": (1, 2), "fs": (3, 1)},
            {"r": 2, "es": (3, 1), "fs": (1, 2)},
            {"r": 2, "es": (3, 2), "fs": (1, 1)},
            {"r": 3, "es": (1, 1, 1), "fs": (3, 1, 1)},
            {"r": 3, "es": (3, 1, 1), "fs": (1, 1, 1)},
            {"r": 3, "es": (1, 1, 1), "fs": (2, 2, 1)},
            {"r": 3, "es": (2, 1, 1), "fs": (1, 2, 1)},
            {"r": 3, "es": (2, 2, 1), "fs": (1, 1, 1)},
            {"r": 4, "es": (1, 1, 1, 1), "fs": (2, 1, 1, 1)},
            {"r": 4, "es": (2, 1, 1, 1), "fs": (1, 1, 1, 1)},
            {"r": 5, "es": (1, 1, 1, 1, 1), "fs": (1, 1, 1, 1, 1)},
        ],
    )
]


@pytest.mark.parametrize("d, exptd", test_cases)
def test_get_isogeny_primes(d, exptd):
    assert splitting_types(d) == exptd
