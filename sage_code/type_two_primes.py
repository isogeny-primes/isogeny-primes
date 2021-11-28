"""type_one_primes.py

    Deals with the Type two primes.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    ====================================================================

"""

import logging


from sage.all import (
    ceil,
    Integer,
    find_root,
    legendre_symbol,
    log,
    pari,
    prime_range,
)  # pylint: disable=no-name-in-module

from .common_utils import x

logger = logging.getLogger(__name__)

GENERIC_UPPER_BOUND = 10 ** 30


def LLS(p):
    return (log(p) + 9 + 2.5 * (log(log(p))) ** 2) ** 2


def get_type_2_uniform_bound(ecdb_type):

    if ecdb_type == "LSS":
        BOUND_TERM = (log(x) + 9 + 2.5 * (log(log(x))) ** 2) ** 2
    elif ecdb_type == "BS":
        # BOUND_TERM = (4*log(x) + 10)**2
        # BOUND_TERM = (3.29*log(x) + 2.96 + 4.9)**2
        BOUND_TERM = (1.881 * log(x) + 2 * 0.34 + 5.5) ** 2
    else:
        raise ValueError("argument must be LSS or BS")

    f = BOUND_TERM ** 6 + BOUND_TERM ** 3 + 1 - x

    try:
        bound = find_root(f, 1000, GENERIC_UPPER_BOUND)
        return ceil(bound)
    except RuntimeError:
        warning_msg = (
            "Warning: Type 2 bound for quadratic field with "
            "discriminant {} failed. Returning generic upper bound"
        ).format(5)
        print(warning_msg)
        return GENERIC_UPPER_BOUND


def get_type_2_bound(K):
    """The bound in the proof of Theorem 6.4 of Larson/Vaintrob, souped up with
    Theorem 5.1 of Bach and Sorenson."""

    # The Bach and Sorenson parameters
    A = 4
    B = 2.5
    C = 5

    n_K = K.degree()
    delta_K = K.discriminant().abs()

    D = 2 * A * n_K
    E = 4 * A * log(delta_K) + 2 * A * n_K * log(12) + 4 * B * n_K + C + 1

    f = x - (D * log(x) + E) ** 4

    try:
        bound = find_root(f, 10, GENERIC_UPPER_BOUND)
        return ceil(bound)
    except RuntimeError:
        warning_msg = (
            "Type 2 bound for quadratic field with "
            "discriminant {} failed. Returning generic upper bound"
        ).format(delta_K)
        logger.warning(warning_msg)
        return GENERIC_UPPER_BOUND


def satisfies_condition_CC(K, p):
    """Determine whether K,p satisfies condition CC.

    Args:
        K ([NumberField]): the number field
        p ([Prime]): the prime p

    Returns: boolean
    """
    for q in prime_range(p / 4):
        for frak_q in K.primes_above(q):
            f = frak_q.residue_class_degree()
            if f % 2 == 1 and q ** f < p / 4:
                if (q ** (2 * f) + q ** f + 1) % p != 0:
                    if legendre_symbol(q, p) == 1:  # i.e. not inert
                        return False
    return True


def satisfies_condition_CC_uniform(possible_odd_f, p):
    """Determine whether degrees,p satisfies condition CC.
    Args:
        K ([NumberField]): the number field
        p ([Prime]): the prime p
    Returns: boolean
    """
    if p % 4 == 1 or p == 2:
        return False
    for q in prime_range((p / 4) ^ (1 / max(possible_odd_f)) + 1):
        if legendre_symbol(q, p) == 1:
            if all((q ** (2 * f) + q ** f + 1) % p != 0 for f in possible_odd_f):
                return False
    return True


def get_type_2_primes(K, bound=None):
    """Compute a list containing the type 2 primes"""

    # First get the bound
    if bound is None:
        bound = get_type_2_bound(K)
        logger.info("type_2_bound = {}".format(bound))

    # We need to include all primes up to 25
    # see Larson/Vaintrob's proof of Theorem 6.4
    output = set(prime_range(25))

    for p in pari.primes(25, bound):
        p_int = Integer(p)
        if p_int % 4 == 3:  # Type 2 primes necessarily congruent to 3 mod 4
            if satisfies_condition_CC(K, p_int):
                output.add(p_int)

    output = list(output)
    output.sort()
    return output
