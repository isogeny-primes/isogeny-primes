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

from sage.all import Partitions, divisors, matrix, GF, lcm
from itertools import product

from .common_utils import get_weil_polys
from .pil_integers import collapse_tuple

def is_non_increasing(t):

    for j in range(len(t)-1):
        if t[j] < t[j+1]:
            return False
    return True


def get_beta_mats_with_pow(F, pow=1):

    frob_polys = get_weil_polys(F)
    return [matrix.companion(f) ** pow for f in frob_polys]


def remove_duplicates(list_of_e_tuples):

    output = []

    for t in list_of_e_tuples:
        if is_non_increasing(t):
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
        for a_part in parts:
            divisors_of_partition = [divisors(t) for t in a_part]
            possible_es = list(product(*divisors_of_partition))
            possible_es = remove_duplicates(possible_es)
            for possible_e_vec in possible_es:
                f_vec = [a_part[j]/possible_e_vec[j] for j in range(r)]
                type_dict = {}
                type_dict['r'] = r
                type_dict['es'] = possible_e_vec
                type_dict['fs'] = tuple(f_vec)
                types.append(type_dict)
    return types


def bound_from_split_type(split_type, eps, q):
    """This implements Algorithm 5.1 from the paper, and returns the integers
       B_{eps,q} and B^*_{eps,q}.

    Args:
        d (int): degree of number field
        eps (tuple): isogeny signature
        q (int): auxiliary rational prime
    """
    frob_poly_mats = [get_beta_mats_with_pow(GF(q**f), pow=my_e) for my_e,f in zip(split_type['es'],split_type['fs'])]
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
            running_lcm = lcm(pil_int, running_lcm)
    if zero_detection_flag:
        return running_lcm, 0
    else:
        return running_lcm, running_lcm


def B_eps_q(d,eps,q):

    split_types = splitting_types(d)
    B_star = 1
    B = 1
    for split_type in split_types:
        pil_int_star, pil_int = bound_from_split_type(split_type, eps, q)
        B_star = lcm(B_star,pil_int_star)
        B = lcm(B,pil_int)
    return B_star, B