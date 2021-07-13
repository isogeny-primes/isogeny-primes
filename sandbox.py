R.<x> = QQ[]
pol = x^4 - x^3 + 4*x^2 + x + 1
K = NumberField(pol,name='a')
Kgal, KtoKgal = K.galois_closure('b', map=True)
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
    C_primes_dict = get_C_primes(Kgal, KtoKgal, G_K, q, epsilons, q_class_group_order, loop_curves=False)
    unified_dict = {}
    q_rat = Integer(q.norm())
    for eps in epsilons:
        unified_dict[eps] = lcm([q_rat, C_primes_dict[eps]])
    tracking_dict[q] = unified_dict

aux_primes_2 = Kgal.primes_above(167)
epsilons_2 = {(12,12,12,12,0,0,0,0) : "blah"}
tracking_dict_2 = {}

# This one is OK
for q in aux_primes_2:
    q_class_group_order = C_Kgal(q).multiplicative_order()
    # these will be dicts with keys the epsilons, values sets of primes
    C_primes_dict = get_C_primes(Kgal, KtoKgal, G_K, q, epsilons_2, q_class_group_order, loop_curves=False)
    unified_dict = {}
    q_rat = Integer(q.norm())
    for eps in epsilons_2:
        unified_dict[eps] = lcm([q_rat, C_primes_dict[eps]])
    tracking_dict_2[q] = unified_dict


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
