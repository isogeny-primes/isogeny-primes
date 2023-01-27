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

def is_non_increasing(t):

    for j in range(len(t)-1):
        if t[j] < t[j+1]:
            return False
    return True


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


def strong_mult_bound(d, eps, q):
    """This implements Algorithm 5.1 from the paper, and returns the integers
       B_{eps,q} and B^*_{eps,q}.

    Args:
        d (int): degree of number field
        eps (tuple): isogeny signature
        q (int): auxiliary rational prime
    """

    frob_polys = get_weil_polys(GF(q))
    frob_poly_mats = [matrix.companion(a_frob_poly) for a_frob_poly in frob_polys]
    matrix_parent = frob_poly_mats[0].parent()
    tr_eps = sum(eps)
    alpha_to_eps_mat = matrix_parent(q**tr_eps)
    B_eps_q = 1
    for a_beta_mat in frob_poly_mats:
        pil_mat = alpha_to_eps_mat.tensor_product(a_beta_mat.parent()(1)) - (
                  alpha_to_eps_mat.parent()(1)).tensor_product(a_beta_mat)
        pil_int = pil_mat.det()
        B_eps_q = lcm(pil_int, B_eps_q)
    return B_eps_q