from itertools import product

from sage.all import prod, ZZ, sqrt, GF

import logging

from .common_utils import x, weil_polynomial_is_elliptic, eps_exp


logger = logging.getLogger(__name__)


def filter_possible_values(possible_values_list, fq, prime_field):

    output = []

    for c in possible_values_list:
        if c ** 2 == prime_field(1):
            output.append(c)
        elif c ** 2 == prime_field(fq.norm()):
            output.append(c)
        else:
            fq_char = ZZ(fq.norm()).prime_divisors()[0]
            possible_mid_coeffs = lifts_in_range(
                2 * sqrt(fq.norm()), c + prime_field(fq.norm()) / c
            )
            possible_weil_polys = [
                x ** 2 + a * x + fq.norm() for a in possible_mid_coeffs
            ]
            # as sanity check, ensure these are Weil polys
            possible_weil_polys = [
                f for f in possible_weil_polys if f.is_weil_polynomial()
            ]

            elliptic_weil_polys = [
                f
                for f in possible_weil_polys
                if weil_polynomial_is_elliptic(f, fq_char, fq.residue_class_degree())
            ]
            if elliptic_weil_polys:
                output.append(c)
    return output


def get_possible_vals_at_gens(gens_info, eps, embeddings, residue_field, prime_field):

    output = {}
    # frak_p0 = K.primes_above(p)[0]  # choice of p_0
    # residue_field = frak_p0.residue_field(names='z')
    # prime_field = GF(p)

    for class_gp_gen in gens_info:
        class_gp_order, alpha = gens_info[class_gp_gen]
        alpha_to_eps = eps_exp(alpha, eps, embeddings)
        alpha_to_eps_mod_p0 = residue_field(alpha_to_eps)
        logger.debug(f"alpha^eps mod p0: {alpha_to_eps_mod_p0}")
        try:
            c_power_12h = prime_field(alpha_to_eps_mod_p0)
        except TypeError:
            # means alpha_to_eps_mod_p0 is not in GF(p) so can ignore p and move on
            return {}
        possible_values = c_power_12h.nth_root(12 * class_gp_order, all=True)
        filtered_values = filter_possible_values(
            possible_values, class_gp_gen, prime_field
        )
        output[class_gp_gen] = list({x ** 12 for x in filtered_values})

    return output


def tuple_exp(tup, exp_tup):
    return tuple((t ** e for t, e in zip(tup, exp_tup)))


def lifts_in_range(N, res_class):

    output = []
    p = res_class.modulus()
    centered_lift = res_class.lift_centered()

    low_run = centered_lift

    while low_run >= -N:
        output.append(low_run)
        low_run = low_run - p

    high_run = centered_lift

    while high_run <= N:
        output.append(high_run)
        high_run = high_run + p

    return [a for a in list(set(output)) if a.abs() <= N]


def get_prime_gens(C_K, my_gens):

    output = []
    it = C_K.number_field().primes_of_bounded_norm_iter(800)

    if not my_gens:
        # means C_K is trivial, so any prime ideal will do
        candidate = next(it)
        output.append(candidate)

    output = [None] * len(my_gens)
    my_gens_copy = set(my_gens)
    while my_gens_copy:
        candidate = next(it)
        candidate_class = C_K(candidate)
        if candidate_class in my_gens:
            output[my_gens.index(candidate_class)] = candidate
            my_gens_copy.remove(candidate_class)

    return output


def final_filter(C_K, Kgal, aux_primes, my_gens_ideals, gens_info, p, eps, embeddings):
    """The possible isogeny prime p is assumed coprime to the prime ideals in my_gens_ideals
    at this point."""
    frak_p0 = Kgal.primes_above(p)[0]  # choice of p_0
    residue_field = frak_p0.residue_field(names="z")
    prime_field = GF(p)

    logger.debug(f"Starting final filter for Prime {p} for eps {eps}")

    # Step 1
    possible_vals_at_gens = get_possible_vals_at_gens(
        gens_info, eps, embeddings, residue_field, prime_field
    )

    if (not possible_vals_at_gens) or (not all(possible_vals_at_gens.values())):
        logger.debug(f"Prime {p} for eps {eps} filtered in Step 1 of Heavy filter")
        logger.debug(f"Possible vals at gens: {possible_vals_at_gens}")
        return False
    logger.debug(f"Possible vals at gens: {possible_vals_at_gens}")
    # Step 2

    # a list of tuples
    possible_vals_cart_prod = list(
        product(*[possible_vals_at_gens[q] for q in my_gens_ideals])
    )

    # Step 3

    # The idea is that we try to filter out each tuple in
    # possible_vals_cart_prod using aux primes; the paper explains how
    # these can be considered as refined epsilon types.

    still_in_the_game = possible_vals_cart_prod.copy()
    for q in aux_primes:
        if q.is_coprime(p):
            exponents_in_class_group = C_K(q).exponents()

            # Check that these exponents correspond to the ideals in
            # my_gens_ideals in the correct order
            sanity_check = prod(
                [Q ** a for Q, a in zip(my_gens_ideals, exponents_in_class_group)]
            )
            assert C_K(sanity_check) == C_K(q)

            the_principal_ideal = q * prod(
                [Q ** (-a) for Q, a in zip(my_gens_ideals, exponents_in_class_group)]
            )
            alphas = the_principal_ideal.gens_reduced()
            assert len(alphas) == 1, "the principal ideal isn't principal!!!"
            alpha = alphas[0]
            alpha_to_eps = eps_exp(alpha, eps, embeddings)
            alpha_to_eps_mod_p0 = residue_field(alpha_to_eps)
            try:
                thingy = prime_field(alpha_to_eps_mod_p0)
            except TypeError as err:
                # means alpha_to_eps_mod_p0 is not in GF(p) so can ignore and move on
                logger.debug(
                    f"Prime {p} for eps {eps} filtered in Step 3a of Heavy filter"
                )
                logger.debug(f"{repr(err)}")
                return False
            new_still_in_the_game = []
            for possible_val in still_in_the_game:
                possible_val_with_raised_exp = tuple_exp(
                    possible_val, exponents_in_class_group
                )
                my_possible_val = [thingy * prod(possible_val_with_raised_exp)]

                filtered_values = filter_possible_values(
                    my_possible_val, q, prime_field
                )
                if filtered_values:
                    new_still_in_the_game.append(possible_val)
            still_in_the_game = new_still_in_the_game
            if not still_in_the_game:
                logger.debug(
                    f"Prime {p} for eps {eps} filtered in Step 3b of Heavy filter"
                )
                return False

    logger.debug(f"Prime {p} for eps {eps} survived Heavy filter")
    # If not returned False by now, then no obstruction to p being an isogeny prime
    return True


def character_enumeration_filter(
    C_K, Kgal, tracking_dict_inv_collapsed, epsilons, aux_primes, embeddings
):
    my_gens = list(C_K.gens())  # done here to fix these, since the order matters later
    my_gens_ideals = get_prime_gens(C_K, my_gens)
    gens_info = {}
    for q in my_gens_ideals:
        q_order = C_K(q).multiplicative_order()
        alphas = (q ** q_order).gens_reduced()
        assert len(alphas) == 1
        alpha = alphas[0]
        gens_info[q] = (q_order, alpha)
    logger.debug(f"Kgal: {Kgal}, C_K: {C_K}")
    logger.debug(f"gen_ideals: {my_gens_ideals}, gen_info: {gens_info}")
    eps_prime_dict = {
        eps: tracking_dict_inv_collapsed[eps].prime_divisors() for eps in epsilons
    }
    all_pre_type_one_two_primes = {p for k in eps_prime_dict for p in eps_prime_dict[k]}
    logger.debug(
        "Pre type one two candidates before filtering: {}".format(
            all_pre_type_one_two_primes
        )
    )
    prime_support_my_gens_ideals = list(
        {a for P in my_gens_ideals for a in ZZ(P.norm()).prime_divisors()}
    )
    eps_prime_filt_dict = {}

    for eps in epsilons:
        survived_primes = []
        for p in eps_prime_dict[eps]:
            if p in prime_support_my_gens_ideals:
                survived_primes.append(p)
            elif final_filter(
                C_K, Kgal, aux_primes, my_gens_ideals, gens_info, p, eps, embeddings
            ):
                survived_primes.append(p)
        eps_prime_filt_dict[eps] = set(survived_primes)

    output = set.union(*(val for val in eps_prime_filt_dict.values()))
    return sorted(output)
