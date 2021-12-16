from itertools import product

from sage.all import GF, ZZ, prod

import logging


from .common_utils import x, weil_polynomial_is_elliptic, eps_exp, primes_iter

logger = logging.getLogger(__name__)


def filter_possible_values(possible_values_list, q, residue_class_degree, prime_field):

    output = []
    fq = q ** residue_class_degree
    for c in possible_values_list:
        if c ** 2 == prime_field(1):
            output.append(c)
        elif c ** 2 == prime_field(fq ** 2):
            output.append(c)
        else:
            possible_mid_coeffs = lifts_in_hasse_range(fq, c + prime_field(fq) / c)
            possible_weil_polys = [x ** 2 + a * x + fq for a in possible_mid_coeffs]

            elliptic_weil_polys = [
                f
                for f in possible_weil_polys
                if weil_polynomial_is_elliptic(f, q, residue_class_degree)
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
        q = class_gp_gen.smallest_integer()
        e = class_gp_gen.residue_class_degree()
        filtered_values = filter_possible_values(possible_values, q, e, prime_field)
        output[class_gp_gen] = list({x ** 12 for x in filtered_values})

    return output


def tuple_exp(tup, exp_tup):
    return tuple((t ** e for t, e in zip(tup, exp_tup)))


def lifts_in_hasse_range(fq, res_class):
    # N = 2 * sqrt(fq)
    # low_run >=-N
    # low_run ** 2 = (-low_run) ** 2 <= N**2 = 4*fq
    fq4 = 4 * fq

    output = []
    p = res_class.modulus()
    centered_lift = res_class.lift_centered()

    low_run = centered_lift

    while low_run ** 2 <= fq4:
        output.append(low_run)
        low_run = low_run - p

    high_run = centered_lift + p

    while high_run ** 2 <= fq4:
        output.append(high_run)
        high_run = high_run + p

    return output


def get_prime_gens(C_K, norm_bound=800):
    gens = list(C_K.gens())
    gen_ideals = list(C_K.gens_ideals())

    it = C_K.number_field().primes_of_bounded_norm_iter(norm_bound)

    prime_gens = [None] * len(gens)
    gens_todo = set(gens)
    while gens_todo:
        candidate = next(it)
        candidate_class = C_K(candidate)
        if candidate_class in gens_todo:
            prime_gens[gens.index(candidate_class)] = candidate
            gens_todo.remove(candidate_class)

    return prime_gens


def character_unit_filter(OK_star_gens, Fp0, eps, embeddings):
    for alpha in OK_star_gens:
        alpha_to_eps = eps_exp(alpha, eps, embeddings)
        alpha_to_eps_mod_p0 = Fp0(alpha_to_eps)
        if alpha_to_eps_mod_p0 != 1:
            return False
    return True


def final_filter(
    C_K,
    Kgal,
    OK_star_gens,
    aux_primes,
    my_gens_ideals,
    gens_info,
    p,
    eps,
    embeddings,
    alpha_cache={},
):
    """The possible isogeny prime p is assumed coprime to the prime ideals in
    my_gens_ideals at this point."""
    frak_p0 = Kgal.prime_above(p)  # choice of p_0
    residue_field = frak_p0.residue_field(names="z")
    prime_field = GF(p)

    logger.debug(f"Starting character enumeration filter for Prime {p} for eps {eps}")

    survived = character_unit_filter(OK_star_gens, residue_field, eps, embeddings)
    if not survived:
        logger.debug(f"Prime {p} for eps {eps} removed using unit filter")
        return survived

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
        # same result as q.is_coprime(p) but faster
        not_is_coprime = p.divides(q.absolute_norm())
        if not_is_coprime:
            continue
        if q in alpha_cache:
            alpha, exponents_in_class_group = alpha_cache[q]
        else:
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
            alpha_cache[q] = (alpha, exponents_in_class_group)
        alpha_to_eps = eps_exp(alpha, eps, embeddings)
        alpha_to_eps_mod_p0 = residue_field(alpha_to_eps)
        try:
            thingy = prime_field(alpha_to_eps_mod_p0)
        except TypeError as err:
            # means alpha_to_eps_mod_p0 is not in GF(p) so can ignore and move on
            logger.debug(f"Prime {p} for eps {eps} filtered in Step 3a of Heavy filter")
            logger.debug(f"{repr(err)}")
            return False
        new_still_in_the_game = []
        for possible_val in still_in_the_game:
            possible_val_with_raised_exp = tuple_exp(
                possible_val, exponents_in_class_group
            )
            my_possible_val = thingy * prod(possible_val_with_raised_exp)
            my_possible_val_roots = my_possible_val.nth_root(12, all=True)
            char = q.smallest_integer()
            e = q.residue_class_degree()
            filtered_values = filter_possible_values(
                my_possible_val_roots, char, e, prime_field
            )
            if filtered_values:
                new_still_in_the_game.append(possible_val)
        still_in_the_game = new_still_in_the_game
        if not still_in_the_game:
            logger.debug(f"Prime {p} for eps {eps} filtered in Step 3b of Heavy filter")
            return False

    logger.debug(f"Prime {p} for eps {eps} survived Heavy filter")
    # If not returned False by now, then no obstruction to p being an isogeny prime
    return True


def character_enumeration_filter(
    K,
    C_K,
    Kgal,
    tracking_dict_inv_collapsed,
    epsilons,
    aux_primes,
    embeddings,
    stop_strategy="auto",
):
    if stop_strategy == "auto":
        enumeration_bound = aux_primes
    OK_star_gens = K.unit_group().gens_values()
    my_gens_ideals = get_prime_gens(C_K)
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
    pre_type_one_two = {p for k in eps_prime_dict for p in eps_prime_dict[k]}
    logger.info(
        f"Pre type one two primes before enumeration filter: {sorted(pre_type_one_two)}"
    )
    prime_support_my_gens_ideals = list(
        {a for P in my_gens_ideals for a in ZZ(P.norm()).prime_divisors()}
    )
    eps_prime_filt_dict = {}

    alpha_cache = {}
    for eps in epsilons:
        survived_primes = []
        for p in eps_prime_dict[eps]:
            if p in prime_support_my_gens_ideals:
                survived_primes.append(p)
                continue
            if stop_strategy == "auto":
                # stop condition:
                # 4sqrt(Nm(q)) > 2p
                # Nm(q) > (p/2)**2
                stop = (p ** 2 // 4) + 1
                if enumeration_bound:
                    stop = min(stop, enumeration_bound)
                aux_primes = primes_iter(K, stop)
            if final_filter(
                C_K,
                Kgal,
                OK_star_gens,
                aux_primes,
                my_gens_ideals,
                gens_info,
                p,
                eps,
                embeddings,
                alpha_cache,
            ):
                survived_primes.append(p)
        eps_prime_filt_dict[eps] = set(survived_primes)

    output = set.union(*(val for val in eps_prime_filt_dict.values()))
    removed = sorted(pre_type_one_two.difference(output))
    logger.info(f"Pre type one two candidates removed by filtering: {removed}")
    return sorted(output)
