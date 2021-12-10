########################################################################
#                                                                      #
#                        PRE TYPE ONE TWO PRIMES                       #
#                                                                      #
########################################################################
import logging
from itertools import product

from sage.all import (
    ZZ,
    GF,
    Integer,
    gcd,
    lcm,
    matrix,
)  # pylint: disable=no-name-in-module
from sage_code.pil_integers import get_PIL_integers

from .character_enumeration import character_enumeration_filter
from .common_utils import (
    R,
    get_weil_polys,
    gal_act_eps,
    eps_exp,
    CLASS_NUMBER_ONE_DISCS,
)

logger = logging.getLogger(__name__)


def get_eps_type(eps):
    """Returns the type of an epsilon (quadratic, quartic, sextic), where
    an epsilon is considered as a tuple
    """

    if 6 in eps:
        if any(t in eps for t in [4, 8]):
            return "mixed"
        return "sextic"
    if any(t in eps for t in [4, 8]):
        if len(set(eps)) == 1:
            # means it's all 4s or all 8s
            return "quartic-diagonal"
        return "quartic-nondiagonal"
    return "quadratic"


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
        eps_orbit = get_redundant_epsilons(
            an_eps, galois_group=galois_group
        )  # dual (and possibly Galois) orbit
        epsilons_output.add(an_eps)
        epsilons.difference_update(eps_orbit)

    return epsilons_output


def get_pre_type_one_two_epsilons(d, galgp=None, heavy_filter=False):
    """This method computes the epsilon group ring characters of Lemma 1 and
    Remark 1 of Momose. The three epsilons of type 1 and 2 are excluded.

    Args:
        d ([int]): Degree of the number field

    Returns:
        dictionary with keys a list of tuples defining the epsilon, and value
        the type of that epsilon
    """

    epsilons_dict = {}

    epsilons_keys = set(product([0, 4, 6, 8, 12], repeat=d))

    epsilons_keys -= {(0,) * d, (6,) * d, (12,) * d}  # remove types 1 and 2 epsilons

    logger.debug("epsilons before filtering: {}".format(len(epsilons_keys)))

    if not heavy_filter:
        epsilons_keys = remove_redundant_epsilons(epsilons_keys, galois_group=galgp)
        logger.debug("epsilons after filtering: {}".format(len(epsilons_keys)))
    else:
        logger.debug("Heavy filtering is on, so no epsilon filtering for now.")

    epsilons_dict = {eps: get_eps_type(eps) for eps in epsilons_keys}

    return epsilons_dict


def _contains_imaginary_quadratic_field_deg_2(K):
    imag_quad = K.is_totally_imaginary()
    hilbert = K.discriminant() in CLASS_NUMBER_ONE_DISCS
    return imag_quad, hilbert


def contains_imaginary_quadratic_field(K):
    """Choosing auxiliary primes in the PreTypeOneTwoCase requires us to
    choose non-principal primes if K contains an imaginary quadratic field."""

    K_deg_abs = K.absolute_degree()

    if K_deg_abs % 2 == 1:
        return False, False

    if K_deg_abs == 2:
        return _contains_imaginary_quadratic_field_deg_2(K)

    quadratic_subfields = K.subfields(2)

    imag_quad_subfields = [
        L for L, _, _ in quadratic_subfields if L.is_totally_imaginary()
    ]

    contains_hilbert_class_field_of_imag_quad = False

    for L in imag_quad_subfields:
        HL = L.hilbert_class_field("c")
        if HL.absolute_degree().divides(K.absolute_degree()):
            K_HL_composite = K.composite_fields(HL)[0]
            if K_HL_composite.absolute_degree() == K_deg_abs:
                contains_hilbert_class_field_of_imag_quad = True
                break

    return (bool(imag_quad_subfields), contains_hilbert_class_field_of_imag_quad)


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


def filter_ABC_primes(K, prime_list, eps_type):
    """Apply congruence and splitting conditions to primes in prime
    list, depending on the type of epsilon

    Args:
        K ([NumberField]): our number field, assumed Galois
        prime_list ([list]): list of primes to filter
        eps_type ([str]): one of 'quadratic', 'quartic', or 'sextic'
    """

    if eps_type == "quadratic":
        # prime must split or ramify in K
        output_list = []

        for p in prime_list:
            if not K.ideal(p).is_prime():
                output_list.append(p)
        return output_list

    if eps_type == "quartic-nondiagonal":
        # prime must split or ramify in K, and be congruent to 2 mod 3
        output_list = []

        for p in prime_list:
            if p % 3 == 2:
                if not K.ideal(p).is_prime():
                    output_list.append(p)
        return output_list

    if eps_type == "quartic-diagonal":
        # prime must be congruent to 2 mod 3
        output_list = []

        for p in prime_list:
            if p % 3 == 2:
                output_list.append(p)
        return output_list

    if eps_type == "sextic":
        # prime must split or ramify in K, and be congruent to 3 mod 4
        output_list = []

        for p in prime_list:
            if p % 4 == 3:
                if not K.ideal(p).is_prime():
                    output_list.append(p)
        return output_list

    if eps_type == "mixed":
        # prime must split or ramify in K, and be congruent to 1 mod 12
        output_list = []

        for p in prime_list:
            if p % 12 == 1:
                if not K.ideal(p).is_prime():
                    output_list.append(p)
        return output_list

    raise ValueError("type must be quadratic, quartic, sextic, or mixed")


def get_aux_primes(K, norm_bound, C_K, h_K, contains_imaginary_quadratic):
    """Get the auxiliary primes, including the emergency aux primes"""

    aux_primes = K.primes_of_bounded_norm(norm_bound)
    completely_split_rat_primes = K.completely_split_primes(B=500)
    if contains_imaginary_quadratic:

        good_primes = [p for p in completely_split_rat_primes if gcd(p, 6 * h_K) == 1]
        list_of_gens = list(C_K.gens())
        i = 0
        while list_of_gens and (i < len(good_primes)):
            a_good_prime = good_primes[i]
            emergency_prime_candidates = K.primes_above(a_good_prime)

            for candidate in emergency_prime_candidates:
                emergency_gen = C_K(candidate)
                if emergency_gen in list_of_gens:
                    if a_good_prime > norm_bound:
                        aux_primes.append(candidate)
                        logger.debug("Emergency aux prime added: {}".format(candidate))
                    list_of_gens.remove(emergency_gen)
            i += 1

        if list_of_gens:
            raise RuntimeError(
                "We have been unable to add enough emergency "
                "auxiliary primes. Try increasing the `B` parameter above."
            )
    else:
        a_good_prime = completely_split_rat_primes[0]
        candidate = K.primes_above(a_good_prime)[0]
        if a_good_prime > norm_bound:
            aux_primes.append(candidate)
            logger.debug("Emergency aux prime added: {}".format(candidate))

    return aux_primes


def get_AB_integers(
    embeddings, frak_q, epsilons, q_class_group_order, multiplicative_bounds={}
):

    output_dict_AB = {}
    alphas = (frak_q ** q_class_group_order).gens_reduced()
    assert len(alphas) == 1, "q^q_class_group_order not principal, which is very bad"
    alpha = alphas[0]
    nm_q = ZZ(frak_q.norm())
    for eps in epsilons:
        alpha_to_eps = eps_exp(alpha, eps, embeddings)
        A = ZZ((alpha_to_eps - 1).norm())
        B = ZZ((alpha_to_eps - (nm_q ** (12 * q_class_group_order))).norm())
        output_dict_AB[eps] = lcm(A, B)
        if multiplicative_bounds.get(eps):
            output_dict_AB[eps] = gcd(output_dict_AB[eps], multiplicative_bounds[eps])
    return output_dict_AB


def alpha_eps_beta_bound(alpha_eps, beta, nm_q_pow_12hq):
    C_mat = alpha_eps ** 2 - alpha_eps * beta.trace() + nm_q_pow_12hq
    N = ZZ(C_mat.det())
    return N


def get_C_integers(
    K,
    embeddings,
    frak_q,
    epsilons,
    q_class_group_order,
    frob_polys,
    multiplicative_bounds={},
):
    nm_q = ZZ(frak_q.absolute_norm())

    # Initialise output dict to empty sets
    output_dict_C = {}
    for eps in epsilons:
        output_dict_C[eps] = 1

    alphas = (frak_q ** q_class_group_order).gens_reduced()
    assert len(alphas) == 1, "q^q_class_group_order not principal, which is very bad"
    alpha = alphas[0]

    alpha_eps_dict = {}
    for eps in epsilons:
        alpha_eps = eps_exp(alpha, eps, embeddings)
        alpha_eps = matrix.companion(alpha_eps.charpoly()).change_ring(ZZ)
        alpha_eps_dict[eps] = alpha_eps

        A = ZZ((alpha_eps - 1).det())
        B = ZZ((alpha_eps - (nm_q ** (12 * q_class_group_order))).det())
        multiplicative_bound = multiplicative_bounds.get(eps)
        if multiplicative_bound:
            A = gcd(A, multiplicative_bound)
            B = gcd(B, multiplicative_bound)
        output_dict_C[eps] = lcm([output_dict_C[eps], A, B])

    nm_q_pow_12hq = nm_q ** (12 * q_class_group_order)
    betas = [
        matrix.companion(frob_poly) ** (12 * q_class_group_order)
        for frob_poly in frob_polys
    ]

    for eps in epsilons:
        alpha_eps = alpha_eps_dict[eps]
        multiplicative_bound = multiplicative_bounds.get(eps)
        for beta in betas:
            N = alpha_eps_beta_bound(alpha_eps, beta, nm_q_pow_12hq)
            if multiplicative_bound:
                N = gcd(N, multiplicative_bound)
            else:
                N = N.abs().perfect_power()[0]
            output_dict_C[eps] = lcm(output_dict_C[eps], N)
            if output_dict_C[eps] == multiplicative_bound:
                break
    return output_dict_C


def get_U_integers(K, epsilons, embeddings):
    """Get divisibilities from the units"""

    unit_gens = K.unit_group().gens_values()
    return {
        eps: gcd(
            [ZZ((eps_exp(u, eps, embeddings) - 1).absolute_norm()) for u in unit_gens]
        )
        for eps in epsilons
    }


def get_pre_type_one_two_primes(
    K, norm_bound=50, loop_curves=False, use_PIL=False, heavy_filter=False
):
    """Pre type 1-2 primes are the finitely many primes outside of which
    the isogeny character is necessarily of type 2 (or 3, which is not relevant
    for us)."""

    (
        contains_imaginary_quadratic,
        contains_hilbert_class_field,
    ) = contains_imaginary_quadratic_field(K)

    if contains_hilbert_class_field:
        raise ValueError(
            "The number field you entered contains the Hilbert "
            "Class field of an imaginary quadratic field. The set "
            "of isogeny primes in this case is therefore infinite."
        )

    # Set up important objects to be used throughout

    Kgal = K.galois_closure("b")
    C_K = K.class_group()
    h_K = C_K.order()
    aux_primes = get_aux_primes(K, norm_bound, C_K, h_K, contains_imaginary_quadratic)
    embeddings = K.embeddings(Kgal)

    # Generate the epsilons

    if K.is_galois():
        G_K = K.galois_group()
        epsilons = get_pre_type_one_two_epsilons(
            K.degree(), galgp=G_K, heavy_filter=heavy_filter
        )
    else:
        epsilons = get_pre_type_one_two_epsilons(K.degree(), heavy_filter=heavy_filter)

    # Now start with the divisibilities. Do the unit computation first

    U_integers_dict = get_U_integers(K, epsilons, embeddings)
    logger.debug("Computed divisibilities from units")

    # Next do the computation of A,B and C integers

    tracking_dict = {}
    frob_polys_dict = {}

    bound_dict = U_integers_dict

    for q in aux_primes:
        q_class_group_order = C_K(q).multiplicative_order()
        nm_q = q.absolute_norm()
        if loop_curves:
            frob_polys = get_weil_polys(GF(nm_q))
        else:
            frob_polys = R.weil_polynomials(2, nm_q)
        frob_polys_dict[q] = frob_polys
        C_integers_dict = get_C_integers(
            Kgal, embeddings, q, epsilons, q_class_group_order, frob_polys, bound_dict
        )
        unified_dict = {}
        q_char = q.smallest_integer()
        for eps in epsilons:
            q_char_eps = gcd(q_char, bound_dict[eps])
            unified_dict[eps] = lcm([q_char_eps, C_integers_dict[eps]])
        bound_dict = unified_dict

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

    # Optionally use the principal ideal lattice for further filtering

    if use_PIL and h_K > 1:
        logger.debug("Using PIL")
        PIL_integers_dict = get_PIL_integers(
            aux_primes, frob_polys_dict, Kgal, epsilons, embeddings, C_K
        )
        for eps in epsilons:
            bound_dict[eps] = ZZ(gcd(bound_dict[eps], PIL_integers_dict[eps]))

    # Split according to epsilon type, get prime divisors, and filter

    if heavy_filter:
        logger.debug("Using Heavy filtering")
        output = character_enumeration_filter(
            K, C_K, Kgal, bound_dict, epsilons, aux_primes, embeddings
        )
        return output

    # Split according to epsilon type, get prime divisors, and filter

    final_split_dict = {}
    for eps_type in set(epsilons.values()):
        eps_type_tracking_dict_inv = {
            eps: ZZ(bound_dict[eps]) for eps in epsilons if epsilons[eps] == eps_type
        }
        eps_type_output = lcm(list(eps_type_tracking_dict_inv.values()))
        if eps_type_output.is_perfect_power():
            eps_type_output = eps_type_output.perfect_power()[0]
        eps_type_output = eps_type_output.prime_divisors()
        eps_type_output = filter_ABC_primes(Kgal, eps_type_output, eps_type)
        final_split_dict[eps_type] = set(eps_type_output)

    # Take union of all primes over all epsilons, sort, and return

    output = set.union(*(val for val in final_split_dict.values()))
    return sorted(output)
