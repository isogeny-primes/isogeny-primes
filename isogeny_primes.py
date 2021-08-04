"""isogeny_primes.py

    Return finite list of isogeny primes attached to a number field.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    ====================================================================

"""

# Imports

import argparse
import json
from pathlib import Path
from itertools import product
import logging
from sage.all import (QQ, next_prime, IntegerRing, prime_range, ZZ, pari,
        PolynomialRing, Integer, Rationals, legendre_symbol, QuadraticField,
        log, exp, find_root, ceil, NumberField, hilbert_class_polynomial,
        RR, EllipticCurve, ModularSymbols, Gamma0, lcm, oo, parent, Matrix,
        gcd, prod, floor, prime_divisors, kronecker_character,
        J0, kronecker_symbol, companion_matrix, euler_phi, DirichletGroup,
        CyclotomicField)


# Global quantitites

GENERIC_UPPER_BOUND = 10**30
EC_Q_ISOGENY_PRIMES = {2,3,5,7,11,13,17,19,37,43,67,163}
CLASS_NUMBER_ONE_DISCS = {-1, -2, -3, -7, -11, -19, -43, -67, -163}
SMALL_GONALITIES = {2,3,5,7,11,13,17,19,23,29,31,37,41,47,59,71}
R = PolynomialRing(Rationals(), 'x')
FORMAL_IMMERSION_DATA_AT_2_PATH = Path('formal_immersion_at_2.json')
BAD_FORMAL_IMMERSION_DATA_PATH = Path('bad_formal_immersion_data.json')
QUADRATIC_POINTS_DATA_PATH = Path('quadratic_points_catalogue.json')

# Global methods

def weil_polynomial_is_elliptic(f,q,a):
    """
    On input of a polynomial f that is a weil polynomial of degree 2 and has constant
    term q^a we check if it actually comes from an elliptic curve over GF(q^a).
    This uses theorem 4.1 of http://archive.numdam.org/article/ASENS_1969_4_2_4_521_0.pdf
    """
    if f[1] % q != 0:
        return True

    if a%2 == 0:
        if f[1] in [-2*q**(a//2), 2*q**(a//2)]:
            return True
        if q%3 != 1 and f[1] in [-q**(a//2), q**(a//2)]:
            return True
        if q%4 != 1 and f[1] == 0:
            return True
    else:
        if q in [2,3]:
            if f[1] in [-q**((a+1)//2),q**((a+1)//2)]:
                return True
        if f[1] == 0:
            return True

    return False

def get_weil_polys(F):
    """
    Returns all degree 2 weil polynomials over F that are actually comming from an elliptic curve.
    """
    q = F.characteristic()
    a = F.degree()
    weil_polys = R.weil_polynomials(2,q**a)
    return [f for f in weil_polys if weil_polynomial_is_elliptic(f,q,a)]


########################################################################
#                                                                      #
#                               WEEDING                                #
#                                                                      #
########################################################################


def oezman_sieve(p,N):
    """Returns True iff p is in S_N. Only makes sense if p ramifies in K"""

    M = QuadraticField(-N)
    h_M = M.class_number()
    H = M.hilbert_class_field('b')
    primes_above_p = M.primes_above(p)

    primes_tot_split_in_hcf = []

    for P in primes_above_p:
        if len(H.primes_above(P)) == h_M:
            primes_tot_split_in_hcf.append(P)

    if not primes_tot_split_in_hcf:
        return False

    f = R(hilbert_class_polynomial(M.discriminant()))
    B = NumberField(f, name='t')
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
    return [chi for chi in H if chi.conductor() == N and chi.multiplicative_order() == K.degree()][0]


def is_torsion_same(p, K, chi, J0_min, B=30, uniform=False):
    """Returns true if the minus part of J0(p) does not gain new torsion when
    base changing to K"""

    d = K.degree()

    if uniform:
        frob_poly_data = [(q, d) for q in prime_range(d+2,B) if q != p]
    else:
        frob_poly_data = [(q, 1) if chi(q) == 1 else (q, d) for q in prime_range(d+2,B) if gcd(q,p) == 1]

    point_counts = []

    for q,i in frob_poly_data:
        frob_pol_q = J0_min.frobenius_polynomial(q)
        frob_mat = companion_matrix(frob_pol_q)
        point_counts.append((frob_mat**i).charpoly()(1))

    # Recall that the rational torsion on J0(p) is entirely contained in
    # the minus part (theorem of Mazur), so checking no-growth of torsion
    # in minus part is done simply as follows

    return J0(p).rational_torsion_order() == gcd(point_counts)


# def is_rank_of_twist_zero(d,S_min):

#     my_map = S_min.rational_period_mapping()
#     tw = M.twisted_winding_element(0,kronecker_character(d))
#     twmap = my_map(tw)
#     return twmap != parent(twmap)(0)

def is_rank_of_twist_zero(chi,ML,S_min_L):
    """Returns true if the rank of the twist of the minus part by the
    character chi is zero"""

    my_map = S_min_L.rational_period_mapping()
    tw = ML.twisted_winding_element(0,chi)
    twmap = my_map(tw)
    return twmap != parent(twmap)(0)


def works_method_of_appendix(p,K):
    """This implements the method of the appendix, returns True if that
    method is able to remove p as an isogeny prime for K."""

    if QuadraticField(-p).class_number() > 2:
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

            if is_torsion_same(p,K,chi,J0_min):
                if is_rank_of_twist_zero(chi,ML,S_min_L):
                    return True
    return False


def apply_quadratic_weeding(candidates, K):
    """Checks whether possible isogeny prime p can be removed for K a
    quadratic field"""

    removed_primes = set()
    Delta_K = K.discriminant()
    D = Delta_K.squarefree_part()
    ramified_primes = Delta_K.prime_divisors()

    with open(QUADRATIC_POINTS_DATA_PATH, 'r') as qdpts_dat_file:
        qdpts_dat = json.load(qdpts_dat_file)

    for p in candidates-EC_Q_ISOGENY_PRIMES:
        if p > 20:
            if str(p) in qdpts_dat:
                data_this_p = qdpts_dat[str(p)]
                if D in data_this_p['known_D']:
                    continue
                if data_this_p['is_complete']:
                    removed_primes.add(p)
                    continue
                removed_p = False
                for q in ramified_primes:
                    if not oezman_sieve(q,p):
                        # Means we have a local obstruction at q
                        logging.debug("Prime {} removed via Oezman sieve".format(p))
                        removed_primes.add(p)
                        removed_p = True
                        break
                if removed_p:
                    continue
                logging.debug("Attempting method of appendix on prime {}".format(p))
                if works_method_of_appendix(p,K):
                    logging.debug("Prime {} removed via method of appendix".format(p))
                    removed_primes.add(p)
            else:
                logging.debug("Attempting method of appendix on prime {}".format(p))
                if works_method_of_appendix(p,K):
                    logging.debug("Prime {} removed via method of appendix".format(p))
                    removed_primes.add(p)
    return removed_primes


def apply_weeding(candidates, K):
    """Wrapper for the methods in this section"""

    if K.degree() == 2:
        return apply_quadratic_weeding(candidates, K)

    elif K.degree().is_prime() and K.is_abelian():
        removed_primes = set()
        for p in candidates-EC_Q_ISOGENY_PRIMES:
            if p > 20:
                logging.debug("Attempting method of appendix on prime {}".format(p))
                if works_method_of_appendix(p,K):
                    logging.debug("Prime {} removed via method of appendix".format(p))
                    removed_primes.add(p)
        return removed_primes

    return set()


########################################################################
#                                                                      #
#                           TYPE ONE PRIMES                            #
#                                                                      #
########################################################################


def R_du(d,u,M,columns=None,a_inv=False):
    """Returns a matrix that can be used to verify formall immersions on X_0(p)
    for all p > 2*M*d, such that p*u = 1 mod M.
    Args:
        d ([int]): degree of number field
        u ([int]): a unit mod M whose formal immersion properties we'd like to check
        M ([int]): an auxilary integer.
    Returns:
        [Matrix]: The Matrix of Corollary 6.8 of Derickx-Kamienny-Stein-Stoll.
    """
    if columns == None:
        columns =[a for a in range(M) if gcd(a,M)==1]
        a_inv=False
    if not a_inv:
        columns = [(a,int((ZZ(1)/a)%M)) for a in columns]
    return Matrix(ZZ,
        [
            [((0 if 2*((r*a[0])%M) < M else 1) -
              (0 if 2*((r*u*a[1])%M) < M else 1)) for a in columns]

        for r in range(1,d+1)])


def get_M(d, M_start=None, M_stop=None, positive_char=True):
    """
    Gets an integer M such that R_du is rank d for all u in (Z/MZ)^*.

    If positive_char=False then R_du only has rank d in characteristic 0
    Otherwise it has rank d in all characteristics > 2
    """
    if not M_start:
        M_start = 3
    if not M_stop:
        # based on trial and error, should be big enough
        # if not we just raise an error
        M_stop=20*d

    for M in range(M_start,M_stop,2):
        columns =[(a,int((ZZ(1)/a)%M)) for a in range(M) if gcd(a,M)==1]
        M_lcm = 1
        for u in range(M):
            if gcd(u,M)!=1:
                continue
            R = R_du(d,u,M,columns,a_inv=True)
            if R.rank() < d:
                break
            assert R.nrows() == d
            elt_divs = R.elementary_divisors()
            if positive_char and elt_divs[-1].prime_to_m_part(2) > 1:
                break
            M_lcm = lcm(M_lcm,elt_divs[-1])
        else:
            return (M, M_lcm)
    raise ValueError("M_stop was to small, no valid value of M < M_stop could be found")


def R_dp(d,p):
    """Return the formal immersion matrix
    Args:
        d ([int]): degree of number field
        p ([prime]): prime whose formal immersion properties we'd like to check
    Returns:
        [Matrix]: The Matrix whose rows are (T_2-3)*T_i e  for i <= d.
                  This is for verifying the linear independance in
                  Corollary 6.4 of Derickx-Kamienny-Stein-Stoll.
    """
    M = ModularSymbols(Gamma0(p),2)
    S = M.cuspidal_subspace()
    S_int = S.integral_structure()
    e = M([0,oo])
    I2 = M.hecke_operator(2)-3
    def get_row(i):
        return S_int.coordinate_vector(S_int(M.coordinate_vector(I2(M.hecke_operator(i)(e)))))
    return Matrix([get_row(i) for i in range(1,d+1)]).change_ring(ZZ)


def is_formall_immersion_fast(d,p):
    """If this function returns true then we have a formall immersion in all characteristics
    > 2. If it returns false then this means nothing.
    """
    R0 = R_du(d,p,2)
    for M in range(3,floor(p/(2*d))):
        u = int((ZZ(1)/p)%M)
        R_M = R_du(d,u,M)
        R0 = R0.augment(R_M)

        divs = R0.elementary_divisors()
        if divs[-1] == 0:
            continue
        if divs[-1].prime_to_m_part(2) == 1:
            return True
    return False

def is_formall_immersion(d,p):
    M = ModularSymbols(Gamma0(p),2)
    S = M.cuspidal_subspace()
    I2 = M.hecke_operator(2)-3
    #assert I2.matrix().rank()==S.dimension()
    D = R_dp(d,p).elementary_divisors()
    if D and D[-1]:
        return int(D[-1].prime_to_m_part(2))
    return 0


def get_bad_formal_immersion_data(d):
    """
    This is the Oesterlé for type 1 primes with modular symbols main routine.
    The computation of get_bad_formal_immersion_data is actually a two step
    rocket. First Proposition 6.8 of Derickx-Kamienny-Stein-Stollis used to
    replace Parents polynomial of degree 6 bound by something reasonable,
    and then Corollary 6.4 is used to go from something reasonable to the exact list.
    """
    assert d > 0

    p_todo = [int(p) for p in prime_range(11)]
    p_done = {}
    q_to_bad_p = {}

    M = get_M(d)[0]

    for p in prime_range(11,2*M*d):
        #first do a relatively cheap test
        if is_formall_immersion_fast(d,p):
            continue
        #this is more expensive
        is_formall = is_formall_immersion(d,p)
        if is_formall:
            if is_formall > 1:
                p_done[int(p)] = is_formall
        else:
            p_todo.append(int(p))

    for p,q_prod in p_done.items():
        for q in prime_divisors(q_prod):
            q_to_bad_p[int(q)] = int(q_to_bad_p.get(q,1)*p)

    return p_todo, q_to_bad_p


def apply_formal_immersion_at_2(output_thus_far, running_prime_dict_2, Kdeg):

    with open(FORMAL_IMMERSION_DATA_AT_2_PATH, 'r') as fi2_dat_file:
        fi2_dat = json.load(fi2_dat_file)

    largest_prime = fi2_dat.pop('largest_prime')

    if not str(Kdeg) in fi2_dat.keys():
        logging.debug("No formal immersion data at 2 with which to filter")
        return output_thus_far

    fi2_this_d = fi2_dat[str(Kdeg)]

    stubborn_set = {p for p in output_thus_far if p < fi2_this_d['smallest_good_formal_immersion_prime']
                                               or p in fi2_this_d['sporadic_bad_formal_immersion_primes']
                                               or p > largest_prime}

    candidate_set = output_thus_far - stubborn_set
    if not candidate_set:
        logging.debug("No candidate primes eligible for formal immersion at 2 filtering")
        return output_thus_far

    prime_divs_at_2 = running_prime_dict_2.prime_divisors()

    output = stubborn_set
    failed_candidates = set()

    for p in candidate_set:
        if p in prime_divs_at_2:
            output.add(p)
        else:
            failed_candidates.add(p)

    logging.debug("Type one primes removed via formal immersion at 2 filtering: {}".format(failed_candidates))

    return output


def get_N(frob_poly, residue_field_card, exponent):
    """Helper method for computing Type 1 primes"""

    if frob_poly.is_irreducible():
        frob_poly_root_field = frob_poly.root_field('a')
    else:
        frob_poly_root_field = IntegerRing()
    roots_of_frob = frob_poly.roots(frob_poly_root_field)
    if len(roots_of_frob) == 1:
        assert roots_of_frob[0][1] == 2
        beta = roots_of_frob[0][0]
        return 1 + residue_field_card ** exponent - 2 * beta ** exponent
    else:
        beta, beta_bar = [r for r,e in roots_of_frob]
        return 1 + residue_field_card ** exponent - beta ** exponent - beta_bar ** exponent


def get_type_1_primes(K, C_K, norm_bound=50, loop_curves=False):
    """Compute the type 1 primes"""

    h_K = C_K.order()

    # Get bad formal immersion data

    if not BAD_FORMAL_IMMERSION_DATA_PATH.is_file():
        logging.debug("No bad formal immersion data found. Computing and adding ...")
        bad_formal_immersion_list, bad_aux_prime_dict = get_bad_formal_immersion_data(K.degree())
        data_for_json_export = {int(K.degree()) : {
                                                "bad_formal_immersion_list" : bad_formal_immersion_list,
                                                "bad_aux_prime_dict" : bad_aux_prime_dict
                                             }
        }
        with open(BAD_FORMAL_IMMERSION_DATA_PATH, 'w') as fp:
            json.dump(data_for_json_export, fp, indent=4)
        logging.debug("Data added")
    else:
        logging.debug("Bad formal immersion data found. Reading to see if it has our data ...")
        with open(BAD_FORMAL_IMMERSION_DATA_PATH, 'r') as bfi_dat_file:
            bfi_dat = json.load(bfi_dat_file)

        if str(K.degree()) in bfi_dat:
            logging.debug("Reading pre-existing data ...")
            bad_formal_immersion_list = bfi_dat[str(K.degree())]['bad_formal_immersion_list']
            bad_aux_prime_dict = bfi_dat[str(K.degree())]['bad_aux_prime_dict']
        else:
            logging.debug("Data not found. Computing new record ...")
            bad_formal_immersion_list, bad_aux_prime_dict = get_bad_formal_immersion_data(K.degree())
            bfi_dat[str(K.degree())] = {
                                        "bad_formal_immersion_list" : bad_formal_immersion_list,
                                        "bad_aux_prime_dict" : bad_aux_prime_dict
                                       }
            with open(BAD_FORMAL_IMMERSION_DATA_PATH, 'w') as fp:
                json.dump(bfi_dat, fp, indent=4)

    aux_primes = prime_range(norm_bound)
    running_prime_dict = {}

    for q in aux_primes:
        frak_q = K.primes_above(q)[0]
        residue_field = frak_q.residue_field(names='z')
        residue_field_card = residue_field.cardinality()
        frak_q_class_group_order = C_K(frak_q).multiplicative_order()
        exponent = 12 * frak_q_class_group_order

        running_primes = q
        if loop_curves:
            weil_polys = get_weil_polys(residue_field)
        else:
            weil_polys = R.weil_polynomials(2, residue_field_card)

        for wp in weil_polys:
            N = get_N(wp, residue_field_card, exponent)
            N = Integer(N)
            if N != 0:
                # else we can ignore since it doesn't arise from an elliptic curve
                running_primes = lcm(running_primes, N)

        if str(q) in bad_aux_prime_dict:
            running_primes = lcm(running_primes, bad_aux_prime_dict[str(q)])

        running_prime_dict[q] = running_primes

    running_prime_dict_2 = running_prime_dict.pop(2)

    output = gcd(list(running_prime_dict.values()))
    output = set(output.prime_divisors())
    output = apply_formal_immersion_at_2(output, running_prime_dict_2, K.degree())
    output = output.union(set(bad_formal_immersion_list))
    Delta_K = K.discriminant().abs()
    output = output.union(set(Delta_K.prime_divisors()))
    third_set = [1+d for d in (12*h_K).divisors()]  # p : (p-1)|12h_K
    output = output.union(set([p for p in third_set if p.is_prime()]))
    output = list(output)
    output.sort()
    return output


########################################################################
#                                                                      #
#                        PRE TYPE ONE TWO PRIMES                       #
#                                                                      #
########################################################################


def group_ring_exp(x, eps, G_K):
    return prod([sigma(x)**my_pow for my_pow, sigma in zip(eps, G_K)])


def gal_act_eps(eps, sigma):
    return tuple(eps[i-1] for i in sigma)


def get_eps_type(eps):
    """Returns the type of an epsilon (quadratic, quartic, sextic), where
    an epsilon is considered as a tuple
    """

    if 6 in eps:
        return 'sextic'
    elif any(t in eps for t in [4,8]):
        return 'quartic'
    else:
        return 'quadratic'


def get_redundant_epsilons(eps, G_action):
    """Redundant epsilons are those in the Galois and dual orbits of a given
    epsilon. They are redundant because they yield the same ABC integers."""

    redundant_epsilons = set()

    for sigma in G_action:
        eps_to_sigma = gal_act_eps(eps, sigma)
        redundant_epsilons.add(eps_to_sigma)
        eps_to_sigma_dual = tuple((12 - x) for x in eps_to_sigma)
        redundant_epsilons.add(eps_to_sigma_dual)

    return redundant_epsilons


def remove_redundant_epsilons(epsilons, G_action):

    epsilons_output = set()

    while epsilons:
        an_eps = epsilons.pop()
        eps_orbit = get_redundant_epsilons(an_eps, G_action)  # galois and dual orbit
        epsilons_output.add(an_eps)
        epsilons.difference_update(eps_orbit)

    return epsilons_output


def get_pre_type_one_two_epsilons(G):
    """This method computes the epsilon group ring characters of Lemma 1 and
    Remark 1 of Momose. The three epsilons of type 1 and 2 are excluded.

    Args:
        d ([int]): Absolute degree of the Galois closure

    Returns:
        dictionary with keys a list of tuples defining the epsilon, and value
        the type of that epsilon
    """

    epsilons_dict = {}
    d = G.order()
    G_action = G.as_finitely_presented_group().as_permutation_group().orbit(tuple(range(1,d+1)), action="OnTuples")

    quartic_epsilons = list(product([0,4,8,12], repeat=d))
    sextic_epsilons = list(product([0,6,12], repeat=d))

    epsilons_keys = set(quartic_epsilons).union(set(sextic_epsilons))

    epsilons_keys -= {(0,)*d, (6,)*d, (12,)*d}  # remove types 1 and 2 epsilons

    logging.debug("epsilons before filtering: {}".format(len(epsilons_keys)))

    epsilons_keys = remove_redundant_epsilons(epsilons_keys, G_action)

    logging.debug("epsilons after filtering: {}".format(len(epsilons_keys)))

    epsilons_dict = {eps: get_eps_type(eps) for eps in epsilons_keys}

    return epsilons_dict


def contains_imaginary_quadratic_field(K):
    """Choosing auxiliary primes in the PreTypeOneTwoCase requires us to
    choose non-principal primes if K contains an imaginary quadratic field."""

    quadratic_subfields = K.subfields(2)

    imag_quad_subfields = [L for L,_,_ in quadratic_subfields if L.is_totally_imaginary()]

    contains_hilbert_class_field_of_imag_quad = False

    for L in imag_quad_subfields:
        HL = L.hilbert_class_field('c')
        if HL.absolute_degree().divides(K.absolute_degree()):
            K_HL_composite = K.composite_fields(HL)[0]
            if K_HL_composite.absolute_degree() == K.absolute_degree():
                contains_hilbert_class_field_of_imag_quad = True
                break

    return (bool(imag_quad_subfields), contains_hilbert_class_field_of_imag_quad)


def filter_ABC_primes(K, prime_list, eps_type):
    """Apply congruence and splitting conditions to primes in prime
    list, depending on the type of epsilon

    Args:
        K ([NumberField]): our number field, assumed Galois
        prime_list ([list]): list of primes to filter
        eps_type ([str]): one of 'quadratic', 'quartic', or 'sextic'
    """

    if eps_type == 'quadratic':
        # no additional restrictions
        return prime_list

    elif eps_type == 'quartic':
        # prime must split or ramify in K, and be congruent to 2 mod 3
        output_list = []

        for p in prime_list:
            if p%3 == 2:
                if not K.ideal(p).is_prime():
                    output_list.append(p)
        return output_list

    elif eps_type == 'sextic':
        # prime must split or ramify in K, and be congruent to 3 mod 4
        output_list = []

        for p in prime_list:
            if p%4 == 3:
                if not K.ideal(p).is_prime():
                    output_list.append(p)
        return output_list

    else:  # should never happen
        raise ValueError("type must be quadratic, quartic, or sextic")


def get_AB_primes(G_K,frak_q,epsilons,q_class_group_order):
    """K is assumed Galois at this point"""

    output_dict_AB = {}
    alphas = (frak_q ** q_class_group_order).gens_reduced()
    assert len(alphas) == 1, "q^q_class_group_order not principal, which is very bad"
    alpha = alphas[0]
    nm_q = ZZ(frak_q.norm())
    for eps in epsilons:
        alpha_to_eps = group_ring_exp(alpha,eps, G_K)
        A = (alpha_to_eps - 1).norm()
        B = (alpha_to_eps - (nm_q ** (12 * q_class_group_order))).norm()
        output_dict_AB[eps] = lcm(A,B)
    return output_dict_AB


def get_C_primes(K, G_K, frak_q, epsilons, q_class_group_order, loop_curves=False):
    """K is assumed Galois at this point, with Galois group G_K"""

    # Initialise output dict to empty sets
    output_dict_C = {}
    for eps in epsilons:
        output_dict_C[eps] = 1

    residue_field = frak_q.residue_field(names='z')
    alphas = (frak_q ** q_class_group_order).gens_reduced()
    assert len(alphas) == 1, "q^q_class_group_order not principal, which is very bad"
    alpha = alphas[0]
    if loop_curves:
        frob_polys_to_loop = get_weil_polys(residue_field)
    else:
        frob_polys_to_loop = R.weil_polynomials(2, residue_field.cardinality())

    for frob_poly in frob_polys_to_loop:
        if frob_poly.is_irreducible():
            frob_poly_root_field = frob_poly.root_field('a')
            _, K_into_KL, L_into_KL, _ = K.composite_fields(frob_poly_root_field, 'c', both_maps=True)[0]
        else:
            frob_poly_root_field = IntegerRing()
        roots_of_frob = frob_poly.roots(frob_poly_root_field)
        betas = [r for r,e in roots_of_frob]

        for beta in betas:
            if beta in K:
                for eps in epsilons:
                    # print('.', end='', flush=True)
                    N = (group_ring_exp(alpha, eps, G_K) - beta ** (12*q_class_group_order)).absolute_norm()
                    N = ZZ(N)
                    output_dict_C[eps] = lcm(output_dict_C[eps], N)
            else:
                for eps in epsilons:
                    # print('.', end='', flush=True)
                    N = (K_into_KL(group_ring_exp(alpha, eps, G_K)) - L_into_KL(beta ** (12*q_class_group_order))).absolute_norm()
                    N = ZZ(N)
                    output_dict_C[eps] = lcm(output_dict_C[eps], N)
    return output_dict_C


def get_pre_type_one_two_primes(K, norm_bound=50, loop_curves=False):
    """Pre type 1-2 primes are the finitely many primes outside of which
    the isogeny character is necessarily of type 2 (or 3, which is not relevant
    for us)."""

    Kgal = K.galois_closure('b')

    contains_imaginary_quadratic, contains_hilbert_class_field = contains_imaginary_quadratic_field(Kgal)

    if contains_hilbert_class_field:
        if K.absolute_degree() < Kgal.absolute_degree():
            # It _would_ be nice to remove the Galois restriction in Momose Lemma 1 ...
            raise NotImplementedError("The Galois closure of the number field you "
                "entered contains the Hilbert Class field of an imaginary quadratic "
                "field. Our algorithm currently cannot deal with this instance.")
        else:
            raise ValueError("The number field you entered contains the Hilbert "
                            "Class field of an imaginary quadratic field. The set "
                            "of isogeny primes in this case is therefore infinite.")

    aux_primes = Kgal.primes_of_bounded_norm(norm_bound)

    try:
        it = Kgal.primes_of_degree_one_iter(max_iterations=1000)
        aux_prime_count = 2
        while aux_prime_count > 0:
            aux_prime_candidate = next(it)
            if (not contains_imaginary_quadratic) or (not aux_prime_candidate.is_principal()):
                if aux_prime_candidate.norm() > norm_bound:
                    aux_primes.append(aux_prime_candidate)
                    logging.info("Emergency aux prime added")
                aux_prime_count -= 1
    except StopIteration:
        # This should never happen, because in this case we are probably containing
        # the Hilbert class field of an imaginary quadratic field and therefore
        # would already have raised above, but I have been unable to prove this.
        raise ValueError("The galois closure of your number field contains an "
                "imaginary quadratic field, and we are unable to find a "
                "non-principal auxiliary prime which is necessary for the algorithm "
                "to work.")

    G_K = Kgal.galois_group()
    C_Kgal = Kgal.class_group()
    logging.debug("C_Kgal = {}".format(C_Kgal))
    epsilons = get_pre_type_one_two_epsilons(G_K)
    tracking_dict = {}

    for q in aux_primes:
        q_class_group_order = C_Kgal(q).multiplicative_order()
        # these will be dicts with keys the epsilons, values sets of primes
        AB_primes_dict = get_AB_primes(G_K,q,epsilons, q_class_group_order)
        C_primes_dict = get_C_primes(Kgal,G_K, q, epsilons, q_class_group_order, loop_curves)
        unified_dict = {}
        q_rat = Integer(q.norm())
        for eps in epsilons:
            unified_dict[eps] = lcm([q_rat, AB_primes_dict[eps], C_primes_dict[eps]])
        tracking_dict[q] = unified_dict

    tracking_dict_inv_collapsed = {}
    for eps in epsilons:
        q_dict = {}
        for q in aux_primes:
            q_dict[q] = tracking_dict[q][eps]
        q_dict_collapsed = gcd(list(q_dict.values()))
        tracking_dict_inv_collapsed[eps] = q_dict_collapsed

    final_split_dict = {}

    for eps_type in set(epsilons.values()):
        eps_type_tracking_dict_inv = {eps:ZZ(tracking_dict_inv_collapsed[eps]) for eps in epsilons if epsilons[eps] == eps_type}
        eps_type_output = lcm(list(eps_type_tracking_dict_inv.values()))
        eps_type_output = eps_type_output.prime_divisors()
        eps_type_output = filter_ABC_primes(Kgal, eps_type_output, eps_type)
        final_split_dict[eps_type] = set(eps_type_output)


    output = set.union(*(val for val in final_split_dict.values()))
    output = list(output)
    output.sort()
    return output


########################################################################
#                                                                      #
#                          TYPE TWO PRIMES                             #
#                                                                      #
########################################################################


def LLS(p):
    return (log(p) + 9 + 2.5 * (log(log(p)))**2)**2


def get_type_2_uniform_bound(ecdb_type):

    if ecdb_type == 'LSS':
        BOUND_TERM = (log(x) + 9 + 2.5 * (log(log(x)))**2)**2
    elif ecdb_type == 'BS':
        # BOUND_TERM = (4*log(x) + 10)**2
        # BOUND_TERM = (3.29*log(x) + 2.96 + 4.9)**2
        BOUND_TERM = (1.881*log(x) + 2*0.34 + 5.5)**2
    else:
        raise ValueError("argument must be LSS or BS")

    f = BOUND_TERM**6 + BOUND_TERM**3 + 1 - x

    try:
        bound = find_root(f,1000,GENERIC_UPPER_BOUND)
        return ceil(bound)
    except RuntimeError:
        warning_msg = ("Warning: Type 2 bound for quadratic field with "
        "discriminant {} failed. Returning generic upper bound").format(5)
        print(warning_msg)
        return GENERIC_UPPER_BOUND


def get_type_2_bound(K):
    """The bound in the proof of Theorem 6.4 of Larson/Vaintrob, souped up with
    Theorem 5.1 of Bach and Sorenson."""

    # The Bach and Sorenson parameters
    A = 4
    B = 2.5
    C = 5

    n_K = K.degree()
    delta_K = K.discriminant().abs()

    D = 2 * A * n_K
    E = 4 * A * log(delta_K) + 2 * A * n_K * log(12) + 4 * B * n_K + C + 1

    x = R.gen()
    f = x - (D*log(x) + E) ** 4

    try:
        bound = find_root(f,10,GENERIC_UPPER_BOUND)
        return ceil(bound)
    except RuntimeError:
        warning_msg = ("Type 2 bound for quadratic field with "
        "discriminant {} failed. Returning generic upper bound").format(delta_K)
        logging.warning(warning_msg)
        return GENERIC_UPPER_BOUND


def satisfies_condition_CC(K,p):
    """Determine whether K,p satisfies condition CC.

    Args:
        K ([NumberField]): the number field
        p ([Prime]): the prime p

    Returns: boolean
    """
    for q in prime_range(p/4):
        for frak_q in K.primes_above(q):
            f = frak_q.residue_class_degree()
            if f%2 == 1 and q**f < p/4:
                if (q**(2*f) + q**f + 1) % p != 0:
                    if legendre_symbol(q,p) == 1:  # i.e. not inert
                        return False
    return True


def satisfies_condition_CC_uniform(possible_odd_f,p):
    """Determine whether degrees,p satisfies condition CC.
    Args:
        K ([NumberField]): the number field
        p ([Prime]): the prime p
    Returns: boolean
    """
    if p%4 == 1 or p==2:
        return False
    for q in prime_range((p/4)^(1/max(possible_odd_f)) + 1):
        if legendre_symbol(q,p) == 1:
            if all((q**(2*f) + q**f + 1) % p != 0 for f in possible_odd_f):
                return False
    return True


def get_type_2_primes(K, bound=None):
    """Compute a list containing the type 2 primes"""

    # First get the bound
    if bound is None:
        bound = get_type_2_bound(K)
        logging.info("type_2_bound = {}".format(bound))

    # We need to include all primes up to 25
    # see Larson/Vaintrob's proof of Theorem 6.4
    output = set(prime_range(25))

    for p in pari.primes(25, bound):
        p_int = Integer(p)
        if p_int % 4 == 3:  # Type 2 primes necessarily congruent to 3 mod 4
            if satisfies_condition_CC(K,p_int):
                output.add(p_int)

    output = list(output)
    output.sort()
    return output


########################################################################
#                                                                      #
#                            DLMV BOUND                                #
#                                                                      #
########################################################################


def DLMV(K):
    """Compute the DLMV bound"""

    # First compute David's C_0

    Delta_K = K.discriminant().abs()
    h_K = K.class_number()
    R_K = K.regulator()
    r_K = K.unit_group().rank()
    delta_K = log(2)/(r_K + 1)
    C_1_K = r_K ** (r_K + 1) * delta_K**(-(r_K - 1)) / 2
    C_2_K = exp(24 * C_1_K * R_K)
    CHEB_DEN_BOUND = (4*log(Delta_K**h_K) + 5*h_K + 5)**2
    C_0 = ((CHEB_DEN_BOUND**(12*h_K))*C_2_K + CHEB_DEN_BOUND**(6*h_K))**4

    # Now the Type 1 and 2 bounds

    type_1_bound = (1 + 3**(12 * h_K))**2
    type_2_bound = get_type_2_bound(K)

    return max(C_0, type_1_bound, type_2_bound)


########################################################################
#                                                                      #
#                      MAIN CALLING FUNCTION                           #
#                                                                      #
########################################################################


def get_isogeny_primes(K, norm_bound, bound=1000, loop_curves=True):

    # Start with some helpful user info

    logging.info("Finding isogeny primes for {}.".format(K))
    logging.info("Bound on auxiliary primes is {}.".format(norm_bound))

    # Get and show PreTypeOneTwoPrimes

    pre_type_one_two_primes = get_pre_type_one_two_primes(K,
                                norm_bound=norm_bound,
                                loop_curves=loop_curves)
    logging.info("pre_type_1_2_primes = {}".format(pre_type_one_two_primes))

    # Get and show TypeOnePrimes

    C_K = K.class_group()

    type_1_primes = get_type_1_primes(K, C_K, norm_bound=norm_bound,
                                         loop_curves=loop_curves)
    logging.info("type_1_primes = {}".format(type_1_primes))

    # Get and show TypeTwoPrimes

    type_2_primes = get_type_2_primes(K, bound=bound)
    logging.info("type_2_primes = {}".format(type_2_primes))

    # Put them all together

    candidates = set.union(set(type_1_primes),
                           set(pre_type_one_two_primes),
                           set(type_2_primes))

    # Try to remove some of these primes via Bruin-Najman and Box tables,
    # Özman sieve, and method of Appendix

    removed_primes = apply_weeding(candidates, K)

    if removed_primes:
        candidates -= removed_primes
        logging.info("Primes removed via weeding = {}".format(removed_primes))
    else:
        logging.debug("No primes removed via weeding")

    return candidates


########################################################################
#                                                                      #
#                            CLI HANDLER                               #
#                                                                      #
########################################################################


def cli_handler(args):

    f = R(args.f)

    K = NumberField(f, name='a')

    loglevel = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S',
                        level=loglevel)
    logging.debug("Debugging level for log messages set.")

    if args.dlmv:
        dlmv_bound = DLMV(K)
        logging.info("DLMV bound for {} is:\n\n{}\n\nwhich is approximately {}".format(K, dlmv_bound, RR(dlmv_bound)))
    else:
        if args.rigorous:
            bound = None
            logging.info("Checking all Type 2 primes up to conjectural bound")
        else:
            bound = args.bound
            logging.warning("Only checking Type 2 primes up to {}. "
                            "To check all, use the PARI/GP script.".format(bound))
        superset = get_isogeny_primes(K, args.norm_bound, bound, args.loop_curves)

        superset_list = list(superset)
        superset_list.sort()
        logging.info("superset = {}".format(superset_list))

        possible_new_isog_primes = superset - EC_Q_ISOGENY_PRIMES
        possible_new_isog_primes_list = list(possible_new_isog_primes)
        possible_new_isog_primes_list.sort()
        logging.info("Possible new isogeny primes = {}".format(possible_new_isog_primes_list))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('f', metavar='f', type=str,
                         help='defining polynomial for the Number field')
    parser.add_argument("--norm_bound", type=int, help="bound on norm of aux primes in PreTypeOneTwo case", default=50)
    parser.add_argument("--loop_curves", action='store_true', help="loop over elliptic curves, don't just loop over all weil polys")
    parser.add_argument("--dlmv", action='store_true', help="get only DLMV bound")
    parser.add_argument("--bound", type=int, help="bound on Type 2 prime search", default=1000)
    parser.add_argument("--rigorous", action='store_true', help="search all Type 2 primes up to conjectural bound")
    parser.add_argument("--verbose", action='store_true', help="get more info printed")

    args = parser.parse_args()
    cli_handler(args)
