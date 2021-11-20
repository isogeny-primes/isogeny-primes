"""weeding.py

    Automated Weeding.

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

# Imports

import json
import logging
from pathlib import Path

from sage.all import (
    J0,
    CyclotomicField,
    DirichletGroup,
    ModularSymbols,
    NumberField,
    QuadraticField,
    companion_matrix,
    euler_phi,
    gcd,
    hilbert_class_polynomial,
    parent,
    prime_range,
)

from .common_utils import EC_Q_ISOGENY_PRIMES, R

logger = logging.getLogger(__name__)

SMALL_GONALITIES = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 59, 71}
QUADRATIC_POINTS_DATA_PATH = Path(
    "sage_code/sdata_files/quadratic_points_catalogue.json"
)


########################################################################
#                                                                      #
#                               WEEDING                                #
#                                                                      #
########################################################################


def oezman_sieve(p, N):
    """Returns True iff p is in S_N. Only makes sense if p ramifies in K"""

    M = QuadraticField(-N)
    h_M = M.class_number()
    H = M.hilbert_class_field("b")
    primes_above_p = M.primes_above(p)

    primes_tot_split_in_hcf = []

    for P in primes_above_p:
        if len(H.primes_above(P)) == h_M:
            primes_tot_split_in_hcf.append(P)

    if not primes_tot_split_in_hcf:
        return False

    f = R(hilbert_class_polynomial(M.discriminant()))
    B = NumberField(f, name="t")
    assert B.degree() == h_M

    possible_nus = B.primes_above(p)

    for nu in possible_nus:
        if nu.residue_class_degree() == 1:
            return True

    return False


def get_dirichlet_character(K):
    """Returns a Dirichlet character whose fixed field is K"""

    N = K.conductor()
    zeta_order = euler_phi(N)  # maybe do this as in LMFDB
    H = DirichletGroup(N, base_ring=CyclotomicField(zeta_order))
    return [
        chi
        for chi in H
        if chi.conductor() == N and chi.multiplicative_order() == K.degree()
    ][0]


def is_torsion_same(p, K, chi, J0_min, B=30, uniform=False):
    """Returns true if the minus part of J0(p) does not gain new torsion when
    base changing to K"""

    d = K.degree()

    if uniform:
        frob_poly_data = [(q, d) for q in prime_range(d + 2, B) if q != p]
    else:
        frob_poly_data = [
            (q, 1) if chi(q) == 1 else (q, d)
            for q in prime_range(d + 2, B)
            if gcd(q, p) == 1
        ]

    point_counts = []

    for q, i in frob_poly_data:
        frob_pol_q = J0_min.frobenius_polynomial(q)
        frob_mat = companion_matrix(frob_pol_q)
        point_counts.append((frob_mat ** i).charpoly()(1))

    # Recall that the rational torsion on J0(p) is entirely contained in
    # the minus part (theorem of Mazur), so checking no-growth of torsion
    # in minus part is done simply as follows

    return J0(p).rational_torsion_order(proof=False) == gcd(point_counts)


def is_rank_of_twist_zero(chi, ML, S_min_L):
    """Returns true if the rank of the twist of the minus part by the
    character chi is zero"""

    my_map = S_min_L.rational_period_mapping()
    tw = ML.twisted_winding_element(0, chi)
    twmap = my_map(tw)
    return twmap != parent(twmap)(0)


def works_method_of_appendix(p, K):
    """This implements the method of the appendix, returns True if that
    method is able to remove p as an isogeny prime for K."""

    if QuadraticField(-p).class_number() > K.degree():
        if p not in SMALL_GONALITIES:
            M = ModularSymbols(p)
            S = M.cuspidal_subspace()
            T = S.atkin_lehner_operator()
            S_min = (T + parent(T)(1)).kernel()
            J0_min = S_min.abelian_variety()

            chi = get_dirichlet_character(K)

            ML = ModularSymbols(p, base_ring=chi.base_ring())
            SL = ML.cuspidal_subspace()
            TL = SL.atkin_lehner_operator()
            S_min_L = (TL + parent(TL)(1)).kernel()

            if is_torsion_same(p, K, chi, J0_min):
                if is_rank_of_twist_zero(chi, ML, S_min_L):
                    return True
    return False


def apply_quadratic_weeding(candidates, K):
    """Checks whether possible isogeny prime p can be removed for K a
    quadratic field"""

    removed_primes = set()
    Delta_K = K.discriminant()
    D = Delta_K.squarefree_part()
    ramified_primes = Delta_K.prime_divisors()

    with open(QUADRATIC_POINTS_DATA_PATH, "r") as qdpts_dat_file:
        qdpts_dat = json.load(qdpts_dat_file)

    for p in candidates - EC_Q_ISOGENY_PRIMES:
        if p > 20:
            if str(p) in qdpts_dat:
                data_this_p = qdpts_dat[str(p)]
                if D in data_this_p["known_D"]:
                    continue
                if data_this_p["is_complete"]:
                    removed_primes.add(p)
                    continue
                removed_p = False
                for q in ramified_primes:
                    if not oezman_sieve(q, p):
                        # Means we have a local obstruction at q
                        logger.debug("Prime {} removed via Oezman sieve".format(p))
                        removed_primes.add(p)
                        removed_p = True
                        break
                if removed_p:
                    continue
                logger.debug("Attempting method of appendix on prime {}".format(p))
                if works_method_of_appendix(p, K):
                    logger.debug("Prime {} removed via method of appendix".format(p))
                    removed_primes.add(p)
            else:
                logger.debug("Attempting method of appendix on prime {}".format(p))
                if works_method_of_appendix(p, K):
                    logger.debug("Prime {} removed via method of appendix".format(p))
                    removed_primes.add(p)
    return removed_primes


def apply_weeding(candidates, K):
    """Wrapper for the methods in this section"""

    if K.degree() == 2:
        return apply_quadratic_weeding(candidates, K)

    if K.degree().is_prime() and K.is_abelian():
        removed_primes = set()
        for p in candidates - EC_Q_ISOGENY_PRIMES:
            if p > 20:
                logger.debug("Attempting method of appendix on prime {}".format(p))
                if works_method_of_appendix(p, K):
                    logger.debug("Prime {} removed via method of appendix".format(p))
                    removed_primes.add(p)
        return removed_primes

    return set()
