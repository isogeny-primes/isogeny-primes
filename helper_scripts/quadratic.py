"""quadratic.py

Runs the main function on several quadratic fields.

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

import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

from sage.all import QuadraticField, Integer
from isogeny_primes import get_isogeny_primes
from sage_code.common_utils import EC_Q_ISOGENY_PRIMES, R, CLASS_NUMBER_ONE_DISCS

R = 50
square_free_D = [D for D in range(-R, R) if Integer(D).is_squarefree() and D != 1]


def main():
    for D in square_free_D:
        K = QuadraticField(D)
        if not K.discriminant() in CLASS_NUMBER_ONE_DISCS:
            superset, _ = get_isogeny_primes(
                K,
                bound=1000,
                ice_filter=True,
                appendix_bound=1000,
                norm_bound=50,
                auto_stop_strategy=True,
                repeat_bound=4,
            )
            possible_new_isog_primes = superset - EC_Q_ISOGENY_PRIMES
            possible_new_isog_primes_list = list(possible_new_isog_primes)
            possible_new_isog_primes_list.sort()
            if possible_new_isog_primes_list:
                print(f"D = {D} possible isogenies = {possible_new_isog_primes_list}")


main()
