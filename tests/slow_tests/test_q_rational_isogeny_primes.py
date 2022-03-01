"""test_q_rational_isogeny_primes.py

Make sure the code runs on Q!

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

from sage.all import QQ, polygen

from isogeny_primes import get_isogeny_primes
from sage_code.common_utils import EC_Q_ISOGENY_PRIMES
import logging


logger = logging.getLogger(__name__)

TEST_SETTINGS = {
    "norm_bound": 20,
    "bound": 1000,
    "ice_filter": True,
}


def test_rational_isogeny_primes():
    x = polygen(QQ)
    K = QQ.extension(x - 1, "D")

    superset, _ = get_isogeny_primes(K, **TEST_SETTINGS)
    improperly_ruled_out = EC_Q_ISOGENY_PRIMES.difference(superset)
    assert improperly_ruled_out == set()
    todo = set(superset).difference(EC_Q_ISOGENY_PRIMES)
    # would be nice if we could automatically do QQ
    # i.e. we could rule out 23 as well.
    assert todo == set([23])
