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

from sage.all import Partitions, divisors, matrix, GF, lcm, gcd, prime_range, prod, Integer
from itertools import product
import logging
from .common_utils import get_weil_polys, x
from .pil_integers import collapse_tuple
from .type_one_primes import cached_bad_formal_immersion_data

logger = logging.getLogger(__name__)

def is_non_increasing(t, part):

    for j in range(len(t)-1):
        if part[j] == part[j+1]:
            if t[j] < t[j+1]:
                return False
    return True


def get_beta_mats_with_pow(F, pow=1):

    frob_polys = get_weil_polys(F)
    output = [matrix.companion(f) ** pow for f in frob_polys]

    # We also need to add \pm 1, \pm q

    f1 = (x - 1) * (x + 1)
    output.append(matrix.companion(f1) ** pow)

    q = F.cardinality()
    fq = (x - q) * (x + q)
    output.append(matrix.companion(fq) ** pow)

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

    for r in range(1,d+1):
        parts = list(Partitions(d, length=r))        
        for part in parts:
            divisors_of_partition = [divisors(t) for t in part]
            possible_es = list(product(*divisors_of_partition))
            possible_es = remove_duplicates(possible_es, part)
            for possible_e_vec in possible_es:
                f_vec = [part[j]/possible_e_vec[j] for j in range(r)]
                type_dict = {}
                type_dict['r'] = r
                type_dict['es'] = possible_e_vec
                type_dict['fs'] = tuple(f_vec)
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
    frob_poly_mats = [get_beta_mats_with_pow(GF(q**f), pow=12*my_e) for my_e,f in zip(split_type['es'],split_type['fs'])]
    beta_mat_tuples = list(product(*frob_poly_mats))
    collapsed_beta_mats = [
        collapse_tuple(a_beta_tuple) for a_beta_tuple in beta_mat_tuples
    ]
    tr_eps = sum(eps)
    q_to_tr_eps = q**tr_eps
    running_lcm = 1
    zero_detection_flag = False
    for a_beta_mat in collapsed_beta_mats:
        matrix_parent = a_beta_mat.parent()
        alpha_to_eps_mat = matrix_parent(q_to_tr_eps)
        pil_mat = alpha_to_eps_mat.tensor_product(a_beta_mat.parent()(1)) - (
                  alpha_to_eps_mat.parent()(1)).tensor_product(a_beta_mat)
        pil_int = pil_mat.det()
        if pil_int == 0:
            zero_detection_flag = True
        else:
            pil_int = gcd(known_mult_bound, pil_int)
            running_lcm = lcm(pil_int, running_lcm)
    if zero_detection_flag:
        return running_lcm, 0
    else:
        return running_lcm, running_lcm


def B_eps_q(d,eps,q, known_mult_bound=0):

    split_types = splitting_types(d)
    B_star = 1
    B = 1
    for split_type in split_types:
        pil_int_star, pil_int = bound_from_split_type(split_type, eps, q, known_mult_bound)
        B_star = gcd(known_mult_bound, lcm(B_star,pil_int_star))
        B = gcd(known_mult_bound, lcm(B,pil_int))
    return B_star, B

def core_loop(d, eps, aux_primes):

    mult_upper_bd = 0
    trace_eps = sum(eps)

    for q in aux_primes:
        B_star, B = B_eps_q(d,eps,q,mult_upper_bd)
        if trace_eps % 6 != 0:
            assert B_star == B
        mult_upper_bd = gcd(mult_upper_bd, B_star)

    return mult_upper_bd


def unif_bd(d, eps):

    tr_eps = sum(eps)

    if tr_eps % 6 != 0:
        aux_primes = prime_range(6)
    elif tr_eps % 12 == 6:
        aux_primes = prime_range(5, 15)
    else:
        raise ValueError("can't deal with this case")

    return core_loop(d, eps, aux_primes)


def type_one_unif_primes(d, q_bd=4):

    assert q_bd > 3

    aux_primes = prime_range(3, q_bd+1)

    mult_upper_bd = 0
    type_one_eps =  (0,) * d
    bad_formal_immersion_list, bad_aux_prime_dict = cached_bad_formal_immersion_data(d)

    for q in aux_primes:
        B_star, B = B_eps_q(d,type_one_eps,q,mult_upper_bd)

        agfi_q = bad_aux_prime_dict.get(str(q),1)

        q_prod = lcm([(q ** f - 1) for f in range(1,d+1)])
        contribution = lcm([B_star, q_prod, agfi_q])
        mult_upper_bd = gcd(mult_upper_bd, contribution)

    output = set(Integer(mult_upper_bd).prime_divisors())
    output = output.union(set(bad_formal_immersion_list))
    return output
    