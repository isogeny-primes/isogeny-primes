import logging
from itertools import product

from sage.all import (
    ZZ,
    gcd,
    lcm,
    matrix,
    prod,
)

from .common_utils import (
    eps_exp,
)

logger = logging.getLogger(__name__)


def collapse_tuple(a_beta_tuple):

    a_beta_tuple_copy = list(a_beta_tuple).copy()
    beta_list = a_beta_tuple_copy
    output = beta_list.pop(0)
    while beta_list:
        output = beta_list.pop(0).tensor_product(output, subdivide=False)
    return output


def as_ZZ_module(G, debug=False):
    """
    Input:
      - An abelian group G

    Output:
      - (H, Z, L) a tripple of ZZ-modules such that:
         1. H is isomorphic to G
         2. Z = ZZ^G.ngens()
         3. L is a submodule of Z such that Z/L=H
         4. The coordinates in H are such that
              G  -> H
              g |-> H(g.exponents())
            is an isomoprhism.
    """
    invs = list(reversed(G.elementary_divisors()))
    if debug:
        assert G.ngens() == len(invs)
        print(invs, [g.order() for g in G.gens()])
        for g, inv in zip(G.gens(), invs):
            assert g.order() == inv
    ZZn = ZZ ** len(invs)
    H = ZZn.submodule(ZZn.gen(i) * invs[i] for i in range(len(invs)))
    return ZZn / H, ZZn, H


def principal_ideal_lattice(aux_primes, class_group, debug=False):
    """
    Input:
      - aux_primes - a list of primes in a numberfield
      - class_group - the classgroup of the same numberfield
    Output:
      - The submodule of ZZ^aux_primes corresponding to the prinicpal ideals
    """
    C_ZZ_mod, C_num, C_den = as_ZZ_module(class_group)
    ZZt = ZZ ** len(aux_primes)
    if debug:
        for q in aux_primes:
            assert prod(g ^ i for g, i in zip(class_group.gens(), class_group(q).exponents())) == class_group(q)
    phi = ZZt.hom(im_gens=[C_num(class_group(q).exponents()) for q in aux_primes], codomain=C_num)
    return phi.inverse_image(C_den)


def get_relevant_beta_mats(aux_primes, relevant_aux_prime_positions, frob_polys_dict):

    output_dict = {}
    for i in relevant_aux_prime_positions:
        do_stuff = [matrix.companion(a_frob_pol) ** 12 for a_frob_pol in frob_polys_dict[aux_primes[i]]]
        output_dict[aux_primes[i]] = do_stuff

    return output_dict


def get_PIL_integers(aux_primes, frob_polys_dict, Kgal, epsilons, embeddings, C_K):

    Lambda = principal_ideal_lattice(aux_primes, C_K)
    Lambda_basis = Lambda.basis()
    logger.debug("Lambda basis = {}".format(Lambda_basis))
    good_basis_elements = [v for v in Lambda_basis if len(v.nonzero_positions()) > 1]
    relevant_aux_prime_positions = {k for v in good_basis_elements for k in v.nonzero_positions()}
    relevant_beta_mats = get_relevant_beta_mats(aux_primes, relevant_aux_prime_positions, frob_polys_dict)

    alphas_dict = {}
    collapsed_beta_mats = {}
    for v in good_basis_elements:
        the_nonzero_positions = v.nonzero_positions()
        alphas = prod([aux_primes[i] ** v[i] for i in the_nonzero_positions]).gens_reduced()
        assert len(alphas) == 1, "uh oh"
        alphas_dict[v] = alphas[0]
        list_list_mats = [relevant_beta_mats[aux_primes[i]] for i in the_nonzero_positions]
        beta_mat_tuples = list(product(*list_list_mats))
        # import pdb; pdb.set_trace()
        collapsed_beta_mats[v] = [collapse_tuple(a_beta_tuple) for a_beta_tuple in beta_mat_tuples]
    logger.debug("Made the alphas and beta_mat_tuples")

    output_dict = {}
    how_many_eps = len(epsilons)
    i = 1
    for eps in epsilons:
        running_gcd = 0
        for v in good_basis_elements:
            running_lcm = 1
            for a_beta_mat in collapsed_beta_mats[v]:
                alpha_to_eps_mat = eps_exp(alphas_dict[v], eps, embeddings).matrix()
                pil_mat = alpha_to_eps_mat.tensor_product(a_beta_mat.parent()(1)) - (
                    alpha_to_eps_mat.parent()(1)
                ).tensor_product(a_beta_mat)
                pil_int = pil_mat.det()
                running_lcm = lcm(pil_int, running_lcm)
            running_gcd = gcd(running_lcm, running_gcd)
        output_dict[eps] = running_gcd
        logger.debug("Successfully computed PIL int for {} epsilons. {} to go".format(i, how_many_eps - i))
        i += 1
    return output_dict
