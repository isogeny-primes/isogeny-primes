"""generic.py

    Deals with the generic signatures.

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

import logging
from itertools import product

from sage.all import (
    ZZ,
    GF,
    gcd,
    lcm,
    matrix,
)  # pylint: disable=no-name-in-module

from .character_enumeration import character_enumeration_filter
from .common_utils import (
    get_weil_polys,
    gal_act_eps,
    eps_exp,
    galois_action_on_embeddings,
    split_primes_iter,
    primes_iter,
    split_embeddings,
    class_group_norm_map,
    get_eps_type,
    filter_ABC_primes,
)

logger = logging.getLogger(__name__)


def get_redundant_epsilons(eps, galois_group=None):
    """Redundant epsilons are those in the dual orbits of a given
    epsilon. They are redundant because they yield the same ABC integers."""

    if galois_group:
        d = galois_group.order()
        G_action = (
            galois_group.as_finitely_presented_group()
            .as_permutation_group()
            .orbit(tuple(range(1, d + 1)), action="OnTuples")
        )

        redundant_epsilons = set()

        for sigma in G_action:
            eps_to_sigma = gal_act_eps(eps, sigma)
            redundant_epsilons.add(eps_to_sigma)
            eps_to_sigma_dual = tuple((12 - x) for x in eps_to_sigma)
            redundant_epsilons.add(eps_to_sigma_dual)
    else:
        redundant_epsilons = {eps, tuple((12 - x) for x in eps)}

    return redundant_epsilons


def remove_redundant_epsilons(epsilons, galois_group=None):

    epsilons_output = set()

    while epsilons:
        an_eps = epsilons.pop()
        eps_orbit = get_redundant_epsilons(an_eps, galois_group=galois_group)  # dual (and possibly Galois) orbit
        epsilons_output.add(an_eps)
        epsilons.difference_update(eps_orbit)

    return epsilons_output


def generic_signatures(d, strong_type_3_epsilons=None, galgp=None, ice_filter=False):
    """Return the generic and weak Type 3 signatures.

    Args:
        d ([int]): Degree of the number field

    Returns:
        dictionary with keys a list of tuples defining the epsilon, and value
        the type of that epsilon
    """

    epsilons_keys = set(product([0, 4, 6, 8, 12], repeat=d))

    epsilons_keys -= {(0,) * d, (6,) * d, (12,) * d}  # remove types 1 and 2 epsilons

    logger.debug("epsilons before filtering: {}".format(len(epsilons_keys)))

    if not ice_filter:
        epsilons_keys = remove_redundant_epsilons(epsilons_keys, galois_group=galgp)
        logger.debug("epsilons after filtering: {}".format(len(epsilons_keys)))
    else:
        logger.debug("ICE filter is on, so no epsilon filtering.")

    if strong_type_3_epsilons is not None:
        actual_type_3_epsilons = set(strong_type_3_epsilons.keys())
        epsilons_keys = epsilons_keys.difference(actual_type_3_epsilons)  # remove strong type 3 epsilons

    epsilons_dict = {eps: get_eps_type(eps) for eps in epsilons_keys}

    return epsilons_dict


def contains_imaginary_quadratic_field(K):
    """Choosing auxiliary primes in the PreTypeOneTwoCase requires us to
    choose non-principal primes if K contains an imaginary quadratic field."""

    K_deg_abs = K.absolute_degree()

    if K_deg_abs % 2 == 1:
        return False

    if K_deg_abs == 2:
        return K.is_totally_imaginary()

    quadratic_subfields = K.subfields(2)

    imag_quad_subfields = [L for L, _, _ in quadratic_subfields if L.is_totally_imaginary()]

    return bool(imag_quad_subfields)


def get_compositum(nf_list, maps=False):

    if len(nf_list) == 1:
        K = nf_list[0]
        return K, [K.hom(K)]

    nf_list_cp = nf_list.copy()

    running_compositum = nf_list_cp.pop(0)

    while nf_list_cp:
        other = nf_list_cp.pop(0)
        running_compositum = running_compositum.composite_fields(other)[0]

    # Now get the maps if requested

    if maps:
        maps_into_compositum = []

        for K in nf_list:
            K_into_comp = K.embeddings(running_compositum)[0]
            maps_into_compositum.append(K_into_comp)

        return running_compositum, maps_into_compositum

    return running_compositum


def get_aux_primes(K, norm_bound, C_K, h_K, contains_imaginary_quadratic):
    """Get the auxiliary primes, including the emergency aux primes"""

    aux_primes = K.primes_of_bounded_norm(norm_bound)
    my_split_prime_iter = split_primes_iter(K)

    a_good_prime = next(my_split_prime_iter)
    candidate = K.primes_above(a_good_prime)[0]
    if a_good_prime > norm_bound:
        aux_primes.append(candidate)
        logger.debug("Emergency aux prime added: {}".format(candidate))

    if contains_imaginary_quadratic and h_K > 1:

        list_of_gens = list(C_K.gens())
        my_primes_iter = primes_iter(K)

        while list_of_gens:

            candidate = next(my_primes_iter)
            candidate_gen = C_K(candidate)
            if candidate_gen in list_of_gens:
                if candidate.norm() > norm_bound:
                    aux_primes.append(candidate)
                    logger.debug(f"Emergency aux prime added: {candidate}")
                list_of_gens.remove(candidate_gen)

    return aux_primes


def alpha_eps_beta_bound(alpha_eps, beta, nm_q_pow_12hq):
    C_mat = alpha_eps**2 - alpha_eps * beta.trace() + nm_q_pow_12hq
    N = ZZ(C_mat.det())
    return N


def ABC_integers(
    embeddings,
    frak_q,
    epsilons,
    q_class_group_order,
    frob_polys,
    multiplicative_bounds=None,
    ensure_C_nonzero=False,
):
    """Computes the ABC integers.

    Args:
        embeddings ([list]): embeddings of K into Kgal
        frak_q ([prime ideal]): the auxiliary prime
        epsilons ([type]): [description]
        q_class_group_order ([Integer]): order of frak_q in C_K
        frob_polys ([list]): Frobenius polynomials to loop over for C
        multiplicative_bounds ([dict], optional): Existing record of
                   multiplicative bounds for each eps. Defaults to None.
        ensure_C_nonzero (bool, optional): Whether to ensure the C integer
         we are building must be nonzero. Used in type_three_not_momose.py.
         Defaults to False.

    Returns:
        [dict]: Dictionary with keys eps, value ABC(eps, frak_q)
    """

    # Some initial setup and precomputation before the main loop

    if multiplicative_bounds is None:
        multiplicative_bounds = {}

    nm_q = ZZ(frak_q.absolute_norm())
    alphas = (frak_q**q_class_group_order).gens_reduced()
    assert len(alphas) == 1, "q^q_class_group_order not principal, which is very bad"
    alpha = alphas[0]
    output_dict = {}
    nm_q_pow_12hq = nm_q ** (12 * q_class_group_order)
    betas = [matrix.companion(frob_poly) ** (12 * q_class_group_order) for frob_poly in frob_polys]

    # The main loop

    for eps in epsilons:

        # Initialise output dict to norm of q

        multiplicative_bound = multiplicative_bounds.get(eps, 0)
        output_dict[eps] = gcd(multiplicative_bound, ZZ(frak_q.smallest_integer()))

        # Compute twisted norm of alpha w.r.t. epsilon

        alpha_eps = eps_exp(alpha, eps, embeddings)
        alpha_eps = matrix.companion(alpha_eps.charpoly()).change_ring(ZZ)

        # Compute A and B integers and add to the LCM after gcd with
        # existing bound

        A = ZZ((alpha_eps - 1).det())
        B = ZZ((alpha_eps - (nm_q ** (12 * q_class_group_order))).det())

        A = gcd(A, multiplicative_bound)
        B = gcd(B, multiplicative_bound)

        output_dict[eps] = lcm([output_dict[eps], A, B])

        # Compute the C integer from Frobenius roots beta

        for beta in betas:
            C = alpha_eps_beta_bound(alpha_eps, beta, nm_q_pow_12hq)
            C = gcd(C, multiplicative_bound)
            C = C.abs().perfect_power()[0]
            if ensure_C_nonzero:
                if C != 0:
                    output_dict[eps] = lcm(output_dict[eps], C)
            else:
                output_dict[eps] = lcm(output_dict[eps], C)
            if output_dict[eps] == multiplicative_bound:
                break
    return output_dict


def get_U_integers(K, epsilons, embeddings):
    """Get divisibilities from the units"""

    unit_gens = K.unit_group().gens_values()
    return {eps: gcd([ZZ((eps_exp(u, eps, embeddings) - 1).absolute_norm()) for u in unit_gens]) for eps in epsilons}


def pre_type_3_class_group_maps(K, embeddings):
    class_group_maps = {}
    for L, phi, _ in K.subfields(degree=2, name="b"):
        if L.discriminant() > 0:
            continue
        split1, _ = split_embeddings(phi, embeddings)
        epsilon1 = tuple(0 if embedding in split1 else 12 for embedding in embeddings)
        epsilon2 = tuple(12 - a for a in epsilon1)
        norm_map = class_group_norm_map(phi, to_C_L=False)
        class_group_maps[epsilon1] = class_group_maps[epsilon2] = norm_map
    return class_group_maps


def get_strong_type_3_epsilons(K, embeddings):

    list_of_gens = list(K.class_group().gens())

    strong_type_3_epsilons = {}
    for L, phi, _ in K.subfields(degree=2, name="b"):
        if L.discriminant() > 0:
            continue
        norm_map = class_group_norm_map(phi, to_C_L=False)
        if all(norm_map(cl) == 0 for cl in list_of_gens):
            # means K contains HCF of L
            split1, _ = split_embeddings(phi, embeddings)
            epsilon1 = tuple(0 if embedding in split1 else 12 for embedding in embeddings)
            epsilon2 = tuple(12 - a for a in epsilon1)
            strong_type_3_epsilons[epsilon1] = L
            strong_type_3_epsilons[epsilon2] = L
    return strong_type_3_epsilons


def generic_primes(
    K,
    norm_bound=50,
    ice_filter=False,
    auto_stop_strategy=True,
    repeat_bound=4,
    character_enumeration_bound=1000,
):
    """Generic primes are the finitely many primes outside of which
    the isogeny character is necessarily of signature type 1, 2 or 3"""

    contains_imaginary_quadratic = contains_imaginary_quadratic_field(K)

    # Set up important objects to be used throughout

    C_K = K.class_group()
    h_K = C_K.order()

    # Generate the epsilons

    if K.is_galois():
        G_K = K.galois_group()
        G_K_emb, _, _, Kgal, embeddings = galois_action_on_embeddings(G_K)
        strong_type_3_epsilons = get_strong_type_3_epsilons(K, embeddings)
        epsilons = generic_signatures(
            K.degree(),
            strong_type_3_epsilons=strong_type_3_epsilons,
            galgp=G_K_emb,
            ice_filter=ice_filter,
        )
    else:
        Kgal = K.galois_closure("b")
        embeddings = K.embeddings(Kgal)
        strong_type_3_epsilons = get_strong_type_3_epsilons(K, embeddings)
        epsilons = generic_signatures(
            K.degree(),
            strong_type_3_epsilons=strong_type_3_epsilons,
            ice_filter=ice_filter,
        )

    # Now start with the divisibilities. Do the unit computation first

    U_integers_dict = get_U_integers(K, epsilons, embeddings)
    logger.debug("Computed divisibilities from units")

    # Next do the computation of A,B and C integers

    frob_polys_dict = {}

    bound_dict = U_integers_dict

    if auto_stop_strategy:
        aux_primes_iter = split_primes_iter(K)
        class_group_maps = pre_type_3_class_group_maps(K, embeddings)
        epsilon_repeats = {eps: repeat_bound for eps in epsilons}
    else:
        aux_primes = get_aux_primes(K, norm_bound, C_K, h_K, contains_imaginary_quadratic)
        aux_primes_iter = aux_primes

    for q in aux_primes_iter:
        epsilons_todo = set(epsilon_repeats) if auto_stop_strategy else epsilons
        q_class_group_order = C_K(q).multiplicative_order()
        nm_q = q.absolute_norm()
        frob_polys = get_weil_polys(GF(nm_q))
        frob_polys_dict[q] = frob_polys
        ABC_integers_dict = ABC_integers(
            embeddings,
            q,
            epsilons_todo,
            q_class_group_order,
            frob_polys,
            bound_dict,
            ensure_C_nonzero=False,
        )
        if auto_stop_strategy:
            for eps in epsilons_todo:
                if ABC_integers_dict[eps] == bound_dict[eps]:
                    if eps in class_group_maps:
                        if class_group_maps[eps](q) != 0:
                            epsilon_repeats[eps] -= 1
                    else:
                        epsilon_repeats[eps] -= 1
                else:
                    epsilon_repeats[eps] = repeat_bound
                if epsilon_repeats[eps] == 0:
                    del epsilon_repeats[eps]

        bound_dict = {**bound_dict, **ABC_integers_dict}
        if auto_stop_strategy and not epsilon_repeats:
            break

    logger.debug(f"Computed generic ABC integers.")
    logger.debug(f"bound dict before enumeration filter: \n{bound_dict}")

    # Barinder you probably don't want me to compute this cause of prime_divisor
    # however I really want to log this from time to time
    """
    import json

    bound_dict_factored = {
        str(eps): str(bound.prime_divisors()) for eps, bound in bound_dict.items()
    }
    logger.info(
        "bound dict before enumeration filter:\n"
        f"{json.dumps(bound_dict_factored, indent=2)}"
    )
    """

    # Filtering stage

    if ice_filter:
        if auto_stop_strategy:
            aux_primes = character_enumeration_bound
        else:
            aux_primes = get_aux_primes(K, norm_bound, C_K, h_K, contains_imaginary_quadratic)
        logger.debug("Using ICE filter")
        output = character_enumeration_filter(
            K,
            C_K,
            Kgal,
            bound_dict,
            epsilons,
            aux_primes,
            embeddings,
            auto_stop_strategy=auto_stop_strategy,
        )

    else:
        final_split_dict = {}
        for eps_type in set(epsilons.values()):
            eps_type_tracking_dict_inv = {eps: ZZ(bound_dict[eps]) for eps in epsilons if epsilons[eps] == eps_type}
            eps_type_output = lcm(list(eps_type_tracking_dict_inv.values()))
            if eps_type_output.is_perfect_power():
                eps_type_output = eps_type_output.perfect_power()[0]
            eps_type_output = eps_type_output.prime_divisors()
            eps_type_output = filter_ABC_primes(K, eps_type_output, eps_type)
            final_split_dict[eps_type] = set(eps_type_output)

        # Take union of all primes over all epsilons
        output = set.union(*(val for val in final_split_dict.values()))

    return strong_type_3_epsilons, embeddings, sorted(output)
