from isogeny_primes import satisfies_condition_CC

def kernel_generator(phi,p1,p2):
    """
    Find a liniar combiantion of p1, p2 such that a*p1+b*p2 is in the kernel of phi
    """
    p = phi.degree()
    assert p.is_prime()
    if phi(p1) == 0:
        return p1
    p0 = p2
    for i in range(p):
        if phi(p0) == 0:
            return p0
        p0+=p1
    raise RuntimeError("Should not reach this point")

def is_type_2(phi):
    """
    On input of an isogeny of prime degree check if it is of type 2
    This is horribly slow, we should probably do this over finite fields.
    """
    p = phi.degree()
    if (p%4) != 3:
        return False
    assert p.is_prime()
    E = phi.domain()
    K = CyclotomicField(p)
    u = GF(p).unit_gens()[0]
    sigma = K.hom([K.0^(u)])
    #sigma is such that theta_p(sigma)=u
    EK = E.change_ring(K)
    p_tors = EK._p_primary_torsion_basis(p,1)
    if len(p_tors) == 0:
        #in this case we can't compute lambda by just looking at the base change to K
        raise NotImplementedError()
    if len(p_tors) == 1:
        P = p_tors[0][0]
    else:
        phi_K = EK.isogeny(phi.kernel_polynomial())
        P = kernel_generator(phi_K,p_tors[0][0],p_tors[1][0])
    lambda_sigma = [i for i in GF(p) if i.lift()*P == P.change_ring(sigma)][0]
    return lambda_sigma^12 == u^6

E11 = [EllipticCurve([1, 1, 1, -305, 7888]),EllipticCurve([0, -1, 1, -887, -10143]),EllipticCurve([1, 1, 0, -3632, 82757])]
phis = [E.isogenies_prime_degree(11)[0] for E in E11]
[[is_type_2(phi),is_type_2(phi.dual()),phi.domain().has_cm()] for phi in phis]

for d in range(1,8):
    CM_discs = sorted({disc for disc,e in cm_orders(d)})
    print(CM_discs)
    for CM_disc in CM_discs:
        print("d, CM_disc", d, CM_disc)
        R.<x> = QQ[]
        f = R(hilbert_class_polynomial(CM_disc))
        M.<j_cm> = QQ.extension(f)

        E = EllipticCurve_from_j(j_cm)

        assert E.has_cm()

        assert E.cm_discriminant() == CM_disc
        for p in (-CM_disc).prime_divisors():
            if p%4 == 1:
                continue
            print(p)
            assert satisfies_condition_CC(M,p)