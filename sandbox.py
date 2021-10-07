R.<x> = QQ[]
pol = x^4 - x^3 + 4*x^2 + x + 1
K = NumberField(pol,name='a')
Kgal, phi = K.galois_closure('b', map=True)
C_K = K.class_group()

aux_primes = Kgal.primes_above(47)
G_K = Kgal.galois_group()
C_Kgal = Kgal.class_group()
epsilons = get_pre_type_one_two_epsilons(G_K)
tracking_dict = {}

# Maarten you probably DON'T want to execute the following for loop, it will take ages and is not necessary
for q in aux_primes:
    q_class_group_order = C_Kgal(q).multiplicative_order()
    # these will be dicts with keys the epsilons, values sets of primes
    C_primes_dict = get_C_primes(Kgal,G_K, q, epsilons, q_class_group_order, loop_curves=False)
    unified_dict = {}
    q_rat = Integer(q.norm())
    for eps in epsilons:
        unified_dict[eps] = lcm([q_rat, C_primes_dict[eps]])
    tracking_dict[q] = unified_dict

aux_primes_2 = K.primes_above(167)
epsilons_2 = {(12,12,12,12,0,0,0,0) : "blah"}
tracking_dict_2 = {}

# This one is OK
for q in aux_primes_2:
    q_class_group_order = C_Kgal(q).multiplicative_order()
    # these will be dicts with keys the epsilons, values sets of primes
    C_primes_dict = get_C_primes(Kgal,embs[0],G_K, q, epsilons_2, q_class_group_order, loop_curves=False)
    unified_dict = {}
    q_rat = Integer(q.norm())
    for eps in epsilons_2:
        unified_dict[eps] = lcm([q_rat, C_primes_dict[eps]])
    tracking_dict_2[q] = unified_dict

# Find a degree one principal ideal of K for which Maarten's code yields 0

def get_C_primes_type_3(q, q_class_group_order):
    """scratch implementation of PreTypeThreePrimes"""
    L, phi, phi_inv = K.subfields(2)[0]
    alphas = (q ** q_class_group_order).gens_reduced()
    assert len(alphas) == 1, "q^q_class_group_order not principal, which is very bad"
    alpha = alphas[0]

    K_over_L = K.relativize(phi,names='g')
    from_rel, to_rel = K_over_L.structure()

    aL = to_rel(alpha).relative_norm()

    residue_field = q.residue_field(names='z')
    frob_polys_to_loop = R.weil_polynomials(2, residue_field.cardinality())

    running_int = 1

    for frob_poly in frob_polys_to_loop:
        if frob_poly.is_irreducible():
            frob_poly_root_field = frob_poly.root_field('a')
            _, L_into_comp, FPRF_into_comp, _ = L.composite_fields(frob_poly_root_field, 'c', both_maps=True)[0]
        else:
            frob_poly_root_field = IntegerRing()
        roots_of_frob = frob_poly.roots(frob_poly_root_field)
        betas = [r for r,e in roots_of_frob]

        for beta in betas:
            if beta in L:
                N = (aL**12 - beta ** (12*q_class_group_order)).absolute_norm()
                N = ZZ(N)
                running_int = lcm(running_int, N)
            else:
                N = (L_into_comp(aL)**12 - FPRF_into_comp(beta) ** (12*q_class_group_order)).absolute_norm()
                N = ZZ(N)
                running_int = lcm(running_int, N)

    return running_int

def get_C_primes(K, phi, G_K, frak_q, epsilons, q_class_group_order, loop_curves=False):
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
            _, K_into_KL, L_into_KL, _ = phi.codomain().composite_fields(frob_poly_root_field, 'c', both_maps=True)[0]
        else:
            frob_poly_root_field = IntegerRing()
        roots_of_frob = frob_poly.roots(frob_poly_root_field)
        betas = [r for r,e in roots_of_frob]

        for beta in betas:
            if beta in K:
                for eps in epsilons:
                    print('.', end='', flush=True)
                    N = (group_ring_exp(phi(alpha), eps, G_K) - beta ** (12*q_class_group_order)).absolute_norm()
                    N = ZZ(N)
                    output_dict_C[eps] = lcm(output_dict_C[eps], N)
            else:
                for eps in epsilons:
                    print('.', end='', flush=True)
                    N = (K_into_KL(group_ring_exp(phi(alpha), eps, G_K)) - L_into_KL(beta ** (12*q_class_group_order))).absolute_norm()
                    N = ZZ(N)
                    output_dict_C[eps] = lcm(output_dict_C[eps], N)
    return output_dict_C

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
                    print('.', end='', flush=True)
                    N = (group_ring_exp(alpha, eps, G_K) - beta ** (12*q_class_group_order)).absolute_norm()
                    N = ZZ(N)
                    output_dict_C[eps] = lcm(output_dict_C[eps], N)
            else:
                for eps in epsilons:
                    print('.', end='', flush=True)
                    N = (K_into_KL(group_ring_exp(alpha, eps, G_K)) - L_into_KL(beta ** (12*q_class_group_order))).absolute_norm()
                    N = ZZ(N)
                    output_dict_C[eps] = lcm(output_dict_C[eps], N)
    return output_dict_C

output = get_C_primes(Kgal, G_K, gal_aux_primes[0], epsilons, 1)

R.<x> = QQ[]
pol = x^4 - x^3 + 4*x^2 + x + 1
K = NumberField(pol,name='a')
Kgal, phi = K.galois_closure('b', map=True)
C_K = K.class_group()
G_K = Kgal.galois_group()
C_Kgal = Kgal.class_group()
epsilons = {(12,12,12,12,0,0,0,0) : "blah"}
embs = K.embeddings(Kgal)

aux_primes = K.primes_above(167)
q = aux_primes[0]

get_C_primes_type_3(q, C_K(q).multiplicative_order())

tracking_dict = {}




# This one is OK
for q in aux_primes:
    q_class_group_order = C_K(q).multiplicative_order()
    # these will be dicts with keys the epsilons, values sets of primes
    C_primes_dict = get_C_primes(Kgal,embs[0],G_K, q, epsilons, q_class_group_order, loop_curves=False)
    unified_dict = {}
    q_rat = Integer(q.norm())
    for eps in epsilons:
        unified_dict[eps] = lcm([q_rat, C_primes_dict[eps]])
    tracking_dict[q] = unified_dict

# This one is OK
K_primes_deg_1_nonprinc = [q for q in K.primes_of_degree_one_list(100) if not q.is_principal()]
for q in K_primes_deg_1_nonprinc:
    q_class_group_order = C_K(q).multiplicative_order()
    # these will be dicts with keys the epsilons, values sets of primes
    a_val = get_C_primes_type_3(q,q_class_group_order)
    assert a_val != 0

len(set(those_vals))



vals = [get_C_primes(Kgal,Phi,G_K, aux_primes[0], epsilons, 1, loop_curves=False)[(12,12,12,12,0,0,0,0)] for Phi in embs]

K_primes_deg_1_princ = [q for q in K.primes_of_degree_one_list(100) if q.is_principal()]

for q in K_primes_deg_1_princ:
    vals = [get_C_primes(Kgal,Phi,G_K, aux_primes[0], epsilons, 1, loop_curves=False)[(12,12,12,12,0,0,0,0)] for Phi in embs]
    vals_set = set(vals)
    if (0 in vals_set) or (len(vals_set) > 1):
        print(q)

import pickle
with open('beta_data.pickle', 'rb') as handle:
    beta_data = pickle.load(handle)

aux_primes = list(beta_data.keys())

x = polygen(QQ);  K.<a> = NumberField(x^3 + 14*x - 1)
aux_primes = K.primes_of_bounded_norm(30)
CK = K.class_group()


R.<x> = QQ[]
f1 = x^3 - 17
f2 = x^4 - 18
f3 = x^5 + 19

K1 = NumberField(f1,'a')
K2 = NumberField(f2,'b')
K3 = NumberField(f3,'c')

my_list = [K1,K2,K3]

the_shtuff = [chr(i+97) for i,F in enumerate(frob_polys_dict[q])]
# the_nfs = [NumberField(F, chr(i+97)) for i,F in enumerate(frob_polys_dict[q])]

x = polygen(QQ);  K.<a> = NumberField(x^3 - 4095)
Kgal = K.galois_closure('b')
C_K = K.class_group()
h_K = C_K.order()
embeddings = K.embeddings(Kgal)
gens_info = {}
R = PolynomialRing(Rationals(), 'x')
my_gens_ideals = C_K.gens_ideals()
for q in my_gens_ideals:
    q_order = C_K(q).multiplicative_order()
    alphas = (q ** q_order).gens_reduced()
    assert len(alphas) == 1
    alpha = alphas[0]
    gens_info[q] = (q_order, alpha)

def eps_exp(x, eps, Sigma):
    return prod([sigma(x)**my_pow for my_pow, sigma in zip(eps, Sigma)])

aux_primes = K.primes_of_bounded_norm(30)
eps = (4,0,6)
p = 17
q = C_K.gens_ideals()[0]
class_gp_order, alpha = gens_info[q]
alpha_to_eps = eps_exp(alpha, eps, embeddings)
frak_p0 = Kgal.primes_above(p)[0]
residue_field = frak_p0.residue_field(names='z')
prime_field = GF(p)
alpha_to_eps_mod_p0 = residue_field(alpha_to_eps)

try:
    c_power_12h = prime_field(alpha_to_eps_mod_p0)
    possible_values = c_power_12h.nth_root(12*class_gp_order, all=True)
    filtered_values = filter_possible_values(possible_values)
except TypeError:
    # means alpha_to_eps_mod_p0 is not in GF(p) so can ignore and move on

def tuple_exp(tup,exp_tup):
    return tuple((t**e for t,e in zip(tup,exp_tup)))

def filter_possible_values(possible_values_list,fq,prime_field):

    output = []

    for c in possible_values_list:
        if c ** 2 == prime_field(1):
            output.append(c)
        elif c ** 2 == prime_field(fq.norm()):
            output.append(c)
        else:
            fq_char = ZZ(fq.norm()).prime_divisors()[0]
            possible_mid_coeffs = [a for a in range(-fq.norm(),fq.norm()+1) if prime_field(a) == c + prime_field(fq.norm()/c)  ]
            possible_weil_polys = [x ** 2 + a * x + fq.norm() for a in possible_mid_coeffs]
            elliptic_weil_polys = [f for f in possible_weil_polys if weil_polynomial_is_elliptic(f,fq_char,fq.residue_class_degree())]
            if elliptic_weil_polys:
                output.append(c)
    return output