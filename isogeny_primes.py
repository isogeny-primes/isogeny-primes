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
from sage.all import (QQ, next_prime, IntegerRing, prime_range, ZZ, pari,
        PolynomialRing, Integer, Rationals, legendre_symbol, QuadraticField,
        log, exp, find_root, ceil, NumberField, hilbert_class_polynomial,
        RR, EllipticCurve, ModularSymbols, Gamma0, lcm, oo, parent, Matrix,
        gcd, prod, floor, prime_divisors)

# Global constants

# The constant which Momose calls Q_2
Q_2 = 3

# Various other Quantities
GENERIC_UPPER_BOUND = 10**30
EC_Q_ISOGENY_PRIMES = {2,3,5,7,11,13,17,19,37,43,67,163}
CLASS_NUMBER_ONE_DISCS = {-1, -2, -3, -7, -11, -19, -43, -67, -163}
R = PolynomialRing(Rationals(), 'x')
FORMAL_IMMERSION_DATA_PATH = Path('bad_formal_immersion_data.json')

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
    Returns als degree 2 weil polynomials over F that are actually comming from an elliptic curve.
    """
    q = F.characteristic()
    a = F.degree()
    weil_polys = R.weil_polynomials(2,q**a)
    return [f for f in weil_polys if weil_polynomial_is_elliptic(f,q,a)]
    


########################################################################
#                                                                      #
#                           ÖZMAN SIEVE                                #
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
        #based on trial and error, should be big enough
        #if not we just raise an error
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
            #print(d,p,u,M,is_R_du_full_rank(d,u,M))
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


def get_type_1_primes(K, C_K, aux_prime_count=3, loop_curves=False):
    """Compute the type 1 primes"""

    h_K = C_K.order()

    # Get bad formal immersion data

    if not FORMAL_IMMERSION_DATA_PATH.is_file():
        print("No bad formal immersion data found. Computing and adding ...")
        bad_formal_immersion_list, bad_aux_prime_dict = get_bad_formal_immersion_data(K.degree())
        data_for_json_export = {int(K.degree()) : {
                                                "bad_formal_immersion_list" : bad_formal_immersion_list,
                                                "bad_aux_prime_dict" : bad_aux_prime_dict
                                             }
        }
        with open(FORMAL_IMMERSION_DATA_PATH, 'w') as fp:
            json.dump(data_for_json_export, fp, indent=4)
        print("Data added")
    else:
        print("Bad formal immersion data found. Reading to see if it has our data ...")
        with open(FORMAL_IMMERSION_DATA_PATH, 'r') as bfi_dat_file:
            bfi_dat = json.load(bfi_dat_file)

        if str(K.degree()) in bfi_dat:
            print("Reading pre-existing data ...")
            bad_formal_immersion_list = bfi_dat[str(K.degree())]['bad_formal_immersion_list']
            bad_aux_prime_dict = bfi_dat[str(K.degree())]['bad_aux_prime_dict']
        else:
            print("Data not found. Computing new record ...")
            bad_formal_immersion_list, bad_aux_prime_dict = get_bad_formal_immersion_data(K.degree())
            bfi_dat[str(K.degree())] = {
                                        "bad_formal_immersion_list" : bad_formal_immersion_list,
                                        "bad_aux_prime_dict" : bad_aux_prime_dict
                                       }
            with open(FORMAL_IMMERSION_DATA_PATH, 'w') as fp:
                json.dump(bfi_dat, fp, indent=4)

    aux_primes = [Q_2]
    prime_to_append = Q_2
    for _ in range(1,aux_prime_count):
        prime_to_append = next_prime(prime_to_append)
        aux_primes.append(prime_to_append)
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

    output = gcd(list(running_prime_dict.values()))
    output = set(output.prime_divisors())
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


def get_pre_type_one_two_epsilons(d):
    """This method computes the epsilon group ring characters of Lemma 1 and
    Remark 1 of Momose. The three epsilons of type 1 and 2 are excluded.

    Args:
        d (int]): Absolute degree of the Galois closure

    Returns:
        dictionary with keys a list of tuples defining the epsilon, and value
        the type of that epsilon
    """

    epsilons_dict = {}

    quartic_epsilons = list(product([0,4,8,12], repeat=d))
    sextic_epsilons = list(product([0,6,12], repeat=d))

    epsilons_keys = set(quartic_epsilons).union(set(sextic_epsilons))

    epsilons_keys -= {(0,)*d, (6,)*d, (12,)*d}  # remove types 1 and 2 epsilons

    epsilons_dict = {eps: get_eps_type(eps) for eps in epsilons_keys}

    return epsilons_dict


def contains_imaginary_quadratic_field(K):
    """Choosing auxiliary primes in the PreTypeOneTwoCase requires us to
    choose non-principal primes if K contains an imaginary quadratic field.
    This method may not be optimal. Maarten, can you do one better?"""

    quadratic_subfields = K.subfields(2)

    for L,_,_ in quadratic_subfields:
        if L.is_totally_imaginary():
            return True

    return False


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


def get_AB_primes(G_K,q,epsilons,q_class_group_order):
    """K is assumed Galois at this point"""

    output_dict_AB = {}
    alphas = (q ** q_class_group_order).gens_reduced()
    assert len(alphas) == 1, "q^q_class_group_order not principal, which is very bad"
    alpha = alphas[0]
    rat_q = ZZ(q.norm())
    assert rat_q.is_prime(), "somehow the degree 1 prime is not prime"
    for eps in epsilons:
        alpha_to_eps = group_ring_exp(alpha,eps, G_K)
        A = (alpha_to_eps - 1).norm()
        B = (alpha_to_eps - (rat_q ** (12 * q_class_group_order))).norm()
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
                    N = (group_ring_exp(alpha, eps, G_K) - beta ** (12*q_class_group_order)).absolute_norm()
                    N = ZZ(N)
                    if N != 0:
                        output_dict_C[eps] = lcm(output_dict_C[eps], N)
            else:
                for eps in epsilons:
                    N = (K_into_KL(group_ring_exp(alpha, eps, G_K)) - L_into_KL(beta ** (12*q_class_group_order))).absolute_norm()
                    N = ZZ(N)
                    if N != 0:
                        output_dict_C[eps] = lcm(output_dict_C[eps], N)
    return output_dict_C


def get_pre_type_one_two_primes(K, aux_prime_count=3, loop_curves=False, verbose_output=False):
    """Pre type 1-2 primes are the finitely many primes outside of which
    the isogeny character is necessarily of type 2 (or 3, which is not relevant
    for us)."""

    tracking_dict = {}
    Kgal = K.galois_closure('b')
    G_K = Kgal.galois_group()
    C_Kgal = Kgal.class_group()
    epsilons = get_pre_type_one_two_epsilons(Kgal.degree())

    if contains_imaginary_quadratic_field(Kgal):
        it = Kgal.primes_of_degree_one_iter()
        aux_primes = []
        while len(aux_primes) < aux_prime_count:
            aux_prime_candidate = next(it)
            if not aux_prime_candidate.is_principal():
                aux_primes.append(aux_prime_candidate)
    else:
        aux_primes = Kgal.primes_of_degree_one_list(aux_prime_count)

    for q in aux_primes:
        q_class_group_order = C_Kgal(q).multiplicative_order()
        # these will be dicts with keys the epsilons, values sets of primes
        AB_primes_dict = get_AB_primes(G_K,q,epsilons, q_class_group_order)
        C_primes_dict = get_C_primes(Kgal,G_K, q, epsilons, q_class_group_order, loop_curves)
        unified_dict = {}
        q_rat = Integer(q.norm())
        assert q_rat.is_prime()
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

    if verbose_output:
        print({eps:ZZ(q_lcm).prime_divisors() for eps,q_lcm in tracking_dict_inv_collapsed.items()})

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
        warning_msg = ("Warning: Type 2 bound for quadratic field with "
        "discriminant {} failed. Returning generic upper bound").format(delta_K)
        print(warning_msg)
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


def get_type_2_primes(K, bound=None):
    """Compute a list containing the type 2 primes"""

    # First get the bound
    if bound is None:
        bound = get_type_2_bound(K)
        print("type_2_bound = {}".format(bound))

    # We need to include all primes up to 25
    # see Larson/Vaintrob's proof of Theorem 6.4
    output = set(prime_range(25))

    for p in pari.primes(25, bound):
        p_int = Integer(p)
        if p_int % 4 == 3:  # Type 2 primes necessarily congruent to 3 mod 4
            if satisfies_condition_CC(K,p_int):
                output.add(p_int)
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


def get_isogeny_primes(K, aux_prime_count, bound=1000, loop_curves=False, verbose_output=False):

    # Start with some helpful user info

    print("\nFinding isogeny primes for {}\n".format(K))
    print("Number of auxiliary primes is {}\n".format(aux_prime_count))

    # Get and show TypeOnePrimes

    C_K = K.class_group()

    type_1_primes = get_type_1_primes(K, C_K, aux_prime_count=aux_prime_count,
                                         loop_curves=loop_curves)
    print("type_1_primes = {}\n".format(type_1_primes))

    # Get and show PreTypeOneTwoPrimes

    pre_type_one_two_primes = get_pre_type_one_two_primes(K,
                                aux_prime_count=aux_prime_count,
                                loop_curves=loop_curves,
                                verbose_output=verbose_output)
    print("pre_type_2_primes = {}\n".format(pre_type_one_two_primes))

    # Get and show TypeTwoPrimes

    type_2_primes = get_type_2_primes(K, bound=bound)
    print("type_2_primes = {}\n".format(type_2_primes))

    # Put them all together and sort the list before returning
    candidates = set.union(set(type_1_primes),
                           set(pre_type_one_two_primes),
                           set(type_2_primes))
    candidates = list(candidates)
    candidates.sort()

    return candidates


########################################################################
#                                                                      #
#                            CLI HANDLER                               #
#                                                                      #
########################################################################


def cli_handler(args):

    f = R(args.f)

    K = NumberField(f, name='a')

    if args.dlmv:
        dlmv_bound = DLMV(K)
        print("DLMV bound for {} is:\n\n{}\n\nwhich is approximately {}".format(K, dlmv_bound, RR(dlmv_bound)))
    else:
        if args.rigorous:
            bound = None
            print("Checking all Type 2 primes up to conjectural bound")
        else:
            bound = args.bound
            print("WARNING: Only checking Type 2 primes up to {}.\n".format(bound))
            print(("To check all, run with '--rigorous', but be advised that "
                "this will take ages and require loads of memory"))
        superset = get_isogeny_primes(K, args.aux_prime_count, bound, args.loop_curves, args.verbose)
        print("superset = {}".format(superset))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('f', metavar='f', type=str,
                         help='defining polynomial for the Number field')
    parser.add_argument("--aux_prime_count", type=int, help="how many auxiliary primes to take", default=5)
    parser.add_argument("--loop_curves", action='store_true', help="loop over elliptic curves, don't just loop over all weil polys")
    parser.add_argument("--dlmv", action='store_true', help="get only DLMV bound")
    parser.add_argument("--bound", type=int, help="bound on Type 2 prime search", default=1000)
    parser.add_argument("--rigorous", action='store_true', help="search all Type 2 primes up to conjectural bound")
    parser.add_argument("--verbose", action='store_true', help="get more info printed")

    args = parser.parse_args()
    cli_handler(args)
