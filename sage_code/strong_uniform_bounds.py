"""strong_uniform_bounds.py

    Computes the strong uniform bounds from the second paper.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2023 Barinder S. Banwait and Maarten Derickx

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

    The authors can be reached at: barinder.s.banwait@gmail.com and
    maarten@mderickx.nl.

    ====================================================================

"""

from sage.all import (
    Partitions,
    divisors,
    matrix,
    GF,
    lcm,
    gcd,
    prime_range,
    prod,
    ZZ,
)
from itertools import product
import logging
from .common_utils import get_weil_polys, x, is_b_smooth
from .pil_integers import collapse_tuple
from .type_one_primes import cached_bad_formal_immersion_data

logger = logging.getLogger(__name__)


def is_non_increasing(t, part):

    for j in range(len(t) - 1):
        if part[j] == part[j + 1]:
            if t[j] < t[j + 1]:
                return False
    return True


def normalize_matrix(m):
    if m.is_scalar():
        m2 = m.submatrix(0, 0, 1, 1)
        m2.set_immutable()
        return m2

    m2 = matrix.companion(m.charpoly()).change_ring(ZZ)
    m2.set_immutable()
    return m2


def get_beta_mats_with_pow(F, pow=1, minimize=True):

    frob_polys = get_weil_polys(F)
    # We also need to add \pm 1, \pm q
    q = F.cardinality()
    frob_polys += [(x - 1), (x + 1), (x - q), (x + q)]
    output = [matrix.companion(f).change_ring(ZZ) ** pow for f in frob_polys]

    if minimize:
        return {normalize_matrix(m) for m in output}

    return output


def remove_duplicates(list_of_e_tuples, part):

    output = []

    for t in list_of_e_tuples:
        if is_non_increasing(t, part):
            output.append(t)
    return output


def splitting_types(d):
    """Implements step 1 of Algorithm 5.1, the possible splitting types

    Args:
        d (int): degree of number field
    """

    types = []

    for r in range(1, d + 1):
        parts = list(Partitions(d, length=r))
        for part in parts:
            divisors_of_partition = [divisors(t) for t in part]
            possible_es = list(product(*divisors_of_partition))
            possible_es = remove_duplicates(possible_es, part)
            for possible_e_vec in possible_es:
                f_vec = [part[j] / possible_e_vec[j] for j in range(r)]
                type_dict = {}
                type_dict["r"] = r
                type_dict["es"] = possible_e_vec
                type_dict["fs"] = tuple(f_vec)
                types.append(type_dict)
    return types


def bound_from_split_type(split_type, eps, q, known_mult_bound=0):
    """This implements Algorithm 5.1 from the paper, and returns the integers
       B_{eps,q} and B^*_{eps,q}.

    Args:
        d (int): degree of number field
        eps (tuple): isogeny signature
        q (int): auxiliary rational prime
    """

    logger.debug(f"Running q={q}, split type={(split_type['es'],split_type['fs'])} for eps={eps}")
    frob_poly_mats = [
        get_beta_mats_with_pow(GF(q**f), pow=12 * my_e) for my_e, f in zip(split_type["es"], split_type["fs"])
    ]
    beta_mat_tuples = list(product(*frob_poly_mats))
    collapsed_beta_mats = [collapse_tuple(a_beta_tuple) for a_beta_tuple in beta_mat_tuples]
    logger.debug(f"len(collapsed_beta_mats)={len(collapsed_beta_mats)}")

    tr_eps = sum(eps)
    q_to_tr_eps = q**tr_eps
    running_lcm = 1
    zero_detection_flag = False
    for a_beta_mat in collapsed_beta_mats:
        matrix_parent = a_beta_mat.parent()

        B_mat = a_beta_mat.parent()(q_to_tr_eps) - a_beta_mat
        B_int = B_mat.det()
        if B_int == 0:

            zero_detection_flag = True
        else:
            B_int = gcd(known_mult_bound, B_int)
            running_lcm = lcm(B_int, running_lcm)
    if zero_detection_flag:
        return running_lcm, 0
    else:
        return running_lcm, running_lcm


def B_eps_q(d, eps, q, known_mult_bound=0):

    split_types = splitting_types(d)
    B_star = 1
    B = 1
    for split_type in split_types:
        pil_int_star, pil_int = bound_from_split_type(split_type, eps, q, known_mult_bound)
        B_star = gcd(known_mult_bound, lcm(B_star, pil_int_star))
        B = gcd(known_mult_bound, lcm(B, pil_int))
    return B_star, B


def core_loop(d, eps, aux_primes):

    mult_upper_bd = 0
    trace_eps = sum(eps)

    for q in aux_primes:
        B_star, B = B_eps_q(d, eps, q, mult_upper_bd)
        if trace_eps % 6 != 0:
            assert B_star == B
        mult_upper_bd = gcd(mult_upper_bd, B_star)
        logger.debug(f"Upperbound after q={q}: {mult_upper_bd}")

    return mult_upper_bd


def unif_bd(d, eps, aux_bound=6):

    tr_eps = sum(eps)

    if tr_eps % 6 == 0:
        raise ValueError(f"can't handle this case since tr eps = {tr_eps} = 0 mod 6")

    aux_primes = prime_range(aux_bound)

    return core_loop(d, eps, aux_primes)


def type_one_unif_bound(d, q_bd=5):

    assert q_bd > 3

    aux_primes = prime_range(3, q_bd + 1)

    mult_upper_bd = 0
    type_one_eps = (0,) * d
    bad_formal_immersion_list, bad_aux_prime_dict = cached_bad_formal_immersion_data(d)

    for q in aux_primes:
        B_star, B = B_eps_q(d, type_one_eps, q, mult_upper_bd)

        agfi_q = bad_aux_prime_dict.get(str(q), 1)

        q_prod = lcm([(q**f - 1) for f in range(1, d + 1)])

        contribution = lcm([B_star, q_prod, agfi_q])
        mult_upper_bd = gcd(mult_upper_bd, contribution)
        logger.debug(f"Upperbound after q={q}: {mult_upper_bd}")

    if logger.isEnabledFor(logging.DEBUG):
        q = 2
        B_star, B = B_eps_q(d, type_one_eps, q, mult_upper_bd)

        q_prod = lcm([(q**f - 1) for f in range(1, d + 1)])
        contribution = lcm([B_star, q_prod])
        mult_upper_bd_2 = gcd(mult_upper_bd, contribution)

        is_smooth, factors = is_b_smooth(mult_upper_bd_2, 10**9)

        factors_str = "{" + ", ".join(str(i) for i in factors) + "}"
        logger.debug(f"Type 1 upperbound if formal immersion works at 2: {factors_str}")

    bad_mult_upper_bd = prod(bad_formal_immersion_list)
    return lcm(mult_upper_bd, bad_mult_upper_bd)
