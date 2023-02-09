"""type_three_not_momose.py

    Deals with the Type 3 signatures of field L in the case where K
    contains the Hilbert class field of L.

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

from sage.all import ZZ, gcd, lcm, GF  # pylint: disable=no-name-in-module

from .common_utils import auxgens, get_weil_polys
from .generic import ABC_integers
from .character_enumeration import character_enumeration_filter


def get_eps_lcm_dict(C_K, epsilons, embeddings, gen_list):

    ABCstar_lcm = {eps: 1 for eps in epsilons}

    for frak_q in gen_list:
        q_class_group_order = C_K(frak_q).multiplicative_order()
        nm_frak_q = frak_q.absolute_norm()
        frob_polys = get_weil_polys(GF(nm_frak_q))
        ABCstar_lcm_one_q = ABC_integers(
            embeddings,
            frak_q,
            epsilons,
            q_class_group_order,
            frob_polys,
            ensure_C_nonzero=True,
        )

        for eps in epsilons:
            ABCstar_lcm[eps] = lcm(ABCstar_lcm[eps], ABCstar_lcm_one_q[eps])

    return ABCstar_lcm


def type_three_not_momose(K, embeddings, strong_type_3_epsilons):
    """Compute a superset of TypeThreeNotMomosePrimes"""

    if len(strong_type_3_epsilons) == 0:
        return [], []

    C_K = K.class_group()
    h_K = C_K.order()

    # Since the data in strong_type_3_epsilons also contains the
    # IQF L, we need to extract only the epsilons for the actual
    # computation. On the other hand, we also want to report the Ls
    # higher up the stack, so we get those as well

    actual_type_3_epsilons = set(strong_type_3_epsilons.keys())
    type_3_fields = set(strong_type_3_epsilons.values())

    if h_K == 1:
        return [], type_3_fields

    aux_gen_list = auxgens(K)

    bound_dict = {eps: 0 for eps in actual_type_3_epsilons}

    for gen_list in aux_gen_list:
        eps_lcm_dict = get_eps_lcm_dict(C_K, actual_type_3_epsilons, embeddings, gen_list)
        for eps in actual_type_3_epsilons:
            bound_dict[eps] = gcd(bound_dict[eps], eps_lcm_dict[eps])

    for eps in actual_type_3_epsilons:
        bound_dict[eps] = lcm(bound_dict[eps], strong_type_3_epsilons[eps].discriminant())

    Kgal = embeddings[0].codomain()
    epsilons_for_ice = {eps: "quadratic-non-constant" for eps in actual_type_3_epsilons}

    output = character_enumeration_filter(
        K,
        C_K,
        Kgal,
        bound_dict,
        epsilons_for_ice,
        1000,
        embeddings,
        auto_stop_strategy=True,
    )

    return output, type_3_fields
