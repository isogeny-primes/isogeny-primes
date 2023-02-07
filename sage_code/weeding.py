"""weeding.py

    Automated Weeding.

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

# Imports

import json
import logging

from sage.all import (
    J0,
    ModularSymbols,
    QuadraticField,
    companion_matrix,
    gcd,
    legendre_symbol,
    oo,
    parent,
    prime_range,
)  # pylint: disable=no-name-in-module

from .common_utils import EC_Q_ISOGENY_PRIMES, SMALL_GONALITIES
from .config import QUADRATIC_POINTS_DATA_PATH

logger = logging.getLogger(__name__)


########################################################################
#                                                                      #
#                               WEEDING                                #
#                                                                      #
########################################################################


def oezman_sieve(p, N):
    """If p is unramified in Q(sqrt(-N)) this always returns True.
    Otherwise returns True iff p is in S_N. Only makes sense if p ramifies in K.
    Also assumes that N is prime (for the time being).
    """
    assert N.is_prime() and N > 2

    M = QuadraticField(-N)
    if p.divides(M.discriminant()):
        return True

    pp = (M * p).factor()[0][0]
    if pp.is_principal():
        if N % 4 == 3:
            return True
        if legendre_symbol(p, N) == 1:
            return True

    return False


def najman_trbovic_filter(unram_primes, ramified_primes):
    """Return True if a possible isogeny prime can be removed via
    Theorem 2.13 of the paper of Najman-Trbovic"""

    absurd_intersection = set(unram_primes).intersection(set(ramified_primes))
    return len(absurd_intersection) > 0


def get_dirichlet_character(K):
    """Returns a Dirichlet character whose fixed field is K"""
    characters = [chi for chi in K.dirichlet_group() if chi.order() == K.degree()]

    if characters:
        return characters[0]
    raise ValueError(f"Galois group of {K} must be cyclic for it to come from a dirichlet character")


def is_torsion_same(p, K, chi, B=30, uniform=False):
    """Returns true if the minus part of J0(p) does not gain new torsion when
    base changing to K"""
    M = ModularSymbols(p)
    S = M.cuspidal_subspace()
    T = S.atkin_lehner_operator()
    S_min = (T + parent(T)(1)).kernel()
    J0_min = S_min.abelian_variety()

    d = K.degree()

    if uniform:
        frob_poly_data = [(q, d) for q in prime_range(d + 2, B) if q != p]
    else:
        frob_poly_data = [(q, 1) if chi(q) == 1 else (q, d) for q in prime_range(d + 2, B) if gcd(q, p) == 1]

    point_counts = []

    for q, i in frob_poly_data:
        frob_pol_q = J0_min.frobenius_polynomial(q)
        frob_mat = companion_matrix(frob_pol_q)
        point_counts.append((frob_mat ** i).charpoly()(1))

    # Recall that the rational torsion on J0(p) is entirely contained in
    # the minus part (theorem of Mazur), so checking no-growth of torsion
    # in minus part is done simply as follows

    return J0(p).rational_torsion_order(proof=False) == gcd(point_counts)


def is_rank_of_twist_zero(p, chi):
    """Returns true if the rank of the twist of the minus part of X_0(p)
    by the character chi is zero"""
    ML = ModularSymbols(p, base_ring=chi.base_ring())
    SL = ML.cuspidal_subspace()
    TL = SL.atkin_lehner_operator()
    S_min_L = (TL + parent(TL)(1)).kernel()

    for S in S_min_L.decomposition():
        my_map = S.rational_period_mapping()
        w = ML([0, oo])
        wmap = my_map(w)
        if wmap == 0:
            return False
        tw = ML.twisted_winding_element(0, chi)
        twmap = my_map(tw)
        if twmap == 0:
            return False

    return True


def works_method_of_appendix(p, K):
    """This implements the method of the appendix, returns True if that
    method is able to remove p as an isogeny prime for K."""

    if QuadraticField(-p).class_number() <= K.degree():
        return False

    if p in SMALL_GONALITIES:
        return False

    chi = get_dirichlet_character(K)

    if K.degree() == 2 and chi(p) == -1:
        return False

    logger.debug("Testing whether torsion is same")
    if is_torsion_same(p, K, chi):
        logger.debug("Torsion is same test passed")
        logger.debug("Testing whether rank is same")
        if is_rank_of_twist_zero(p, chi):
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
                unram_primes = data_this_p["unramified_primes"]
                if najman_trbovic_filter(unram_primes, ramified_primes):
                    logger.debug("Prime {} removed via Najman-Trbovic filter".format(p))
                    removed_primes.add(p)
    return removed_primes


def apply_weeding(candidates, K, appendix_bound=1000):
    """Wrapper for the methods in this section"""

    removed_primes = set()

    if K.degree() == 2:
        removed_primes = apply_quadratic_weeding(candidates, K)

    if K.degree().is_prime() and K.is_abelian():

        if K.conductor() > 1e6:

            logger.warning(
                "The conductor of your field is too large "
                "for method of appendix to work in reasonable time, "
                "so not doing it."
            )
            return removed_primes

        for p in candidates - EC_Q_ISOGENY_PRIMES - removed_primes:
            if 20 < p < appendix_bound:
                logger.debug(f"Attempting method of appendix on prime {p}")
                if works_method_of_appendix(p, K):
                    logger.debug(f"Prime {p} removed via method of appendix")
                    removed_primes.add(p)
                else:
                    logger.debug(f"Method of appendix fails for {p}")

    return removed_primes
