"""lmfdb_cubic.py

Connects to the LMFDB and runs the main function on some cubic fields
satisfying some parameters.

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

import requests
from sage.all import NumberField
from isogeny_primes import get_isogeny_primes
from sage_code.common_utils import EC_Q_ISOGENY_PRIMES, R

URL_TRUNK = (
    "https://www.lmfdb.org/api/nf_fields/?_format=json&degree=3&"
    "is_galois={}&class_number={}&r2=0&_fields=label,coeffs"
)


def main(is_galois, class_number):
    the_url = URL_TRUNK.format(str(is_galois).lower(), str(class_number))
    payload = requests.get(url=the_url).json()["data"]
    for dat in payload[:10]:
        poly = R(dat["coeffs"])
        label = dat["label"]
        print(label)

        K = NumberField(poly, name="a")
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
        msg = "poly = {} ; label = {} ; possible isogenies = {}"
        print(msg.format(poly, label, possible_new_isog_primes_list))


main(False, 1)
