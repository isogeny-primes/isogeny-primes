"""type_one_primes.py

    Deals with the Type one primes.

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

import json
import logging
from typing import Set

from sage.all import (
    ZZ,
    Gamma0,
    Matrix,
    ModularSymbols,
    floor,
    gcd,
    lcm,
    oo,
    prime_divisors,
    prime_range,
)
from sage.rings.finite_rings.finite_field_constructor import GF

from .common_utils import R, get_weil_polys
from .config import FORMAL_IMMERSION_DATA_AT_2_PATH, BAD_FORMAL_IMMERSION_DATA_PATH

logger = logging.getLogger(__name__)


########################################################################
#                                                                      #
#                           TYPE ONE PRIMES                            #
#                                                                      #
########################################################################


def R_du(d, u, M, columns=None, a_inv=False):
    """Returns a matrix that can be used to verify formall immersions on X_0(p)
    for all p > 2*M*d, such that p*u = 1 mod M.
    Args:
        d ([int]): degree of number field
        u ([int]): a unit mod M whose formal immersion properties we'd like to check
        M ([int]): an auxilary integer.
    Returns:
        [Matrix]: The Matrix of Corollary 6.8 of Derickx-Kamienny-Stein-Stoll.
    """
    if columns is None:
        columns = [a for a in range(M) if gcd(a, M) == 1]
        a_inv = False
    if not a_inv:
        columns = [(a, int((ZZ(1) / a) % M)) for a in columns]
    return Matrix(
        ZZ,
        [
            [
                (
                    (0 if 2 * ((r * a[0]) % M) < M else 1)
                    - (0 if 2 * ((r * u * a[1]) % M) < M else 1)
                )
                for a in columns
            ]
            for r in range(1, d + 1)
        ],
    )


def get_M(d, M_start=None, M_stop=None, positive_char=True):
    """
    Gets an integer M such that R_du is rank d for all u in (Z/MZ)^*.

    If positive_char=False then R_du only has rank d in characteristic 0
    Otherwise it has rank d in all characteristics > 2.
    """
    if not M_start:
        M_start = 3
    if not M_stop:
        # based on trial and error, should be big enough
        # if not we just raise an error
        M_stop = 20 * d

    for M in range(M_start, M_stop, 2):
        columns = [(a, int((ZZ(1) / a) % M)) for a in range(M) if gcd(a, M) == 1]
        M_lcm = 1
        for u in range(M):
            if gcd(u, M) != 1:
                continue
            r_du = R_du(d, u, M, columns, a_inv=True)
            if r_du.rank() < d:
                break
            assert r_du.nrows() == d
            elt_divs = r_du.elementary_divisors()
            if positive_char and elt_divs[-1].prime_to_m_part(2) > 1:
                break
            M_lcm = lcm(M_lcm, elt_divs[-1])
        else:
            return (M, M_lcm)
    raise ValueError("M_stop was to small, no valid value of M < M_stop could be found")


def R_dp(d, p):
    """Return the formal immersion matrix
    Args:
        d ([int]): degree of number field
        p ([prime]): prime whose formal immersion properties we'd like to check
    Returns:
        [Matrix]: The Matrix whose rows are (T_2-3)*T_i e  for i <= d.
                  This is for verifying the linear independance in
                  Corollary 6.4 of Derickx-Kamienny-Stein-Stoll.
    """
    M = ModularSymbols(Gamma0(p), 2)
    S = M.cuspidal_subspace()
    S_int = S.integral_structure()
    e = M([0, oo])
    I2 = M.hecke_operator(2) - 3

    def get_row(i):
        return S_int.coordinate_vector(
            S_int(M.coordinate_vector(I2(M.hecke_operator(i)(e))))
        )

    return Matrix([get_row(i) for i in range(1, d + 1)]).change_ring(ZZ)


def is_formall_immersion_fast(d, p):
    """If this function returns true then we have a formall immersion in all
    characteristics > 2. If it returns false then this means nothing.
    """
    R0 = R_du(d, p, 2)
    for M in range(3, floor(p / (2 * d))):
        u = int((ZZ(1) / p) % M)
        R_M = R_du(d, u, M)
        R0 = R0.augment(R_M)

        divs = R0.elementary_divisors()
        if divs[-1] == 0:
            continue
        if divs[-1].prime_to_m_part(2) == 1:
            return True
    return False


def is_formall_immersion(d, p):
    D = R_dp(d, p).elementary_divisors()
    if D and D[-1]:
        return int(D[-1].prime_to_m_part(2))
    return 0


def get_bad_formal_immersion_data(d):
    """
    This is the Oesterlé for type 1 primes with modular symbols main routine.
    The computation is actually a two step rocket. First Proposition 6.8 of
    Derickx-Kamienny-Stein-Stoll is used to replace Parents polynomial of
    degree 6 bound by something reasonable, and then Corollary 6.4 is used
    to go from something reasonable to the exact list.
    """
    assert d > 0

    p_todo = [int(p) for p in prime_range(11)]
    p_done = {}
    q_to_bad_p = {}

    M = get_M(d)[0]

    for p in prime_range(11, 2 * M * d):
        # first do a relatively cheap test
        if is_formall_immersion_fast(d, p):
            continue
        # this is more expensive
        is_formall = is_formall_immersion(d, p)
        if is_formall:
            if is_formall > 1:
                p_done[int(p)] = is_formall
        else:
            p_todo.append(int(p))

    for p, q_prod in p_done.items():
        for q in prime_divisors(q_prod):
            q_to_bad_p[int(q)] = int(q_to_bad_p.get(q, 1) * p)

    return p_todo, q_to_bad_p


def apply_formal_immersion_at_2(
    output_thus_far: Set[int], running_prime_dict_2: int, Kdeg: int
):

    with open(FORMAL_IMMERSION_DATA_AT_2_PATH, "r") as fi2_dat_file:
        fi2_dat = json.load(fi2_dat_file)

    largest_prime = fi2_dat.pop("largest_prime")

    if not str(Kdeg) in fi2_dat.keys():
        logger.debug("No formal immersion data at 2 with which to filter")
        return output_thus_far

    fi2_this_d = fi2_dat[str(Kdeg)]

    stubborn_set = {
        p
        for p in output_thus_far
        if p < fi2_this_d["smallest_good_formal_immersion_prime"]
        or p in fi2_this_d["sporadic_bad_formal_immersion_primes"]
        or p > largest_prime
    }

    candidate_set = output_thus_far - stubborn_set
    if not candidate_set:
        logger.debug("No candidate primes eligible for formal immersion at 2 filtering")
        return output_thus_far

    output = stubborn_set
    failed_candidates = set()

    for p in candidate_set:
        if p.divides(running_prime_dict_2):
            output.add(p)
        else:
            failed_candidates.add(p)

    logger.debug(
        "Type one primes removed via formal immersion at 2 filtering: {}".format(
            failed_candidates
        )
    )

    return output


def get_N(frob_poly, nm_q, exponent):
    """Helper method for computing Type 1 primes"""
    beta = Matrix.companion(frob_poly) ** exponent
    N = ZZ(1 - beta.trace() + nm_q ** exponent)
    return N


def get_C_integer_type1(K, q, bad_aux_prime_dict, C_K, bound_so_far, loop_curves):
    running_primes = gcd(q, bound_so_far)
    if str(q) in bad_aux_prime_dict:
        running_primes = lcm(
            running_primes, gcd(bad_aux_prime_dict[str(q)], bound_so_far)
        )

    norms_clexp = {
        (frak_q.absolute_norm(), C_K(frak_q).multiplicative_order())
        for frak_q in K.primes_above(q)
    }
    for nm_q, frak_q_class_group_order in norms_clexp:

        exponent = 12 * frak_q_class_group_order

        N_cusp = ZZ(nm_q) ** exponent - 1
        N_cusp = gcd(N_cusp, bound_so_far)
        running_primes = lcm(running_primes, N_cusp)

        if loop_curves:
            weil_polys = get_weil_polys(GF(nm_q))
        else:
            weil_polys = R.weil_polynomials(2, nm_q)

        for wp in weil_polys:
            N = get_N(wp, nm_q, exponent)
            assert N != 0
            N = gcd(N, bound_so_far)

            running_primes = lcm(running_primes, N)
            if running_primes == bound_so_far:
                break

    return gcd(running_primes, bound_so_far)


def get_type_1_primes(K, C_K, norm_bound=50, loop_curves=False):
    """Compute the type 1 primes"""

    h_K = C_K.order()

    # Get bad formal immersion data

    if not BAD_FORMAL_IMMERSION_DATA_PATH.is_file():
        logger.debug("No bad formal immersion data found. Computing and adding ...")
        bad_formal_immersion_list, bad_aux_prime_dict = get_bad_formal_immersion_data(
            K.degree()
        )
        data_for_json_export = {
            int(K.degree()): {
                "bad_formal_immersion_list": bad_formal_immersion_list,
                "bad_aux_prime_dict": bad_aux_prime_dict,
            }
        }
        with open(BAD_FORMAL_IMMERSION_DATA_PATH, "w") as fp:
            json.dump(data_for_json_export, fp, indent=4)
        logger.debug("Data added")
    else:
        logger.debug(
            "Bad formal immersion data found. Reading to see if it has our data ..."
        )
        with open(BAD_FORMAL_IMMERSION_DATA_PATH, "r") as bfi_dat_file:
            bfi_dat = json.load(bfi_dat_file)

        if str(K.degree()) in bfi_dat:
            logger.debug("Reading pre-existing data ...")
            bad_formal_immersion_list = bfi_dat[str(K.degree())][
                "bad_formal_immersion_list"
            ]
            bad_aux_prime_dict = bfi_dat[str(K.degree())]["bad_aux_prime_dict"]
        else:
            logger.debug("Data not found. Computing new record ...")
            (
                bad_formal_immersion_list,
                bad_aux_prime_dict,
            ) = get_bad_formal_immersion_data(K.degree())
            bfi_dat[str(K.degree())] = {
                "bad_formal_immersion_list": bad_formal_immersion_list,
                "bad_aux_prime_dict": bad_aux_prime_dict,
            }
            with open(BAD_FORMAL_IMMERSION_DATA_PATH, "w") as fp:
                json.dump(bfi_dat, fp, indent=4)

    aux_primes = prime_range(3, norm_bound + 1)
    bound_so_far = 0

    for q in aux_primes:
        bound_so_far = get_C_integer_type1(
            K, q, bad_aux_prime_dict, C_K, bound_so_far, loop_curves
        )

    bound_at_2 = get_C_integer_type1(
        K, 2, bad_aux_prime_dict, C_K, bound_so_far, loop_curves
    )

    output = set(bound_so_far.prime_divisors())
    output = apply_formal_immersion_at_2(output, bound_at_2, K.degree())
    output = output.union(set(bad_formal_immersion_list))

    return sorted(output)
