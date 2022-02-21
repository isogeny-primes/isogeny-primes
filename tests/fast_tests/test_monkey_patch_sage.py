"""test_monkey_patch_sage.py

Test the monkey patch.

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

from sage.misc.misc_c import prod
from sage.rings.number_field.number_field import QuadraticField

from sage_code.monkey_patch_sage import ideal_log_relation


def test_ideal_log_relation():
    K = QuadraticField(-41)
    q = (K * 3).factor()[0][0]
    C_K = K.class_group()
    t = ideal_log_relation(q)
    e = q.ideal_class_log()
    assert q == t * prod(gi ** ei.sage() for gi, ei in zip(C_K.gens_ideals(), e))
