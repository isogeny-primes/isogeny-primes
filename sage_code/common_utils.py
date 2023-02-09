"""common_utils.py

    Common methods and classes for multiple parts of the routine.

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

from sage.all import PolynomialRing, Rationals, prod, oo
from sage.arith.misc import primes
from sage.combinat.permutation import Permutation
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.perm_gps.permgroup_named import TransitiveGroup, SymmetricGroup

R = PolynomialRing(Rationals(), "x")
x = R.gen()
EC_Q_ISOGENY_PRIMES = {2, 3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163}
CLASS_NUMBER_ONE_DISCS = {-3, -4, -7, -8, -11, -19, -43, -67, -163}
SMALL_GONALITIES = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 59, 71}

# Global methods


def weil_polynomial_is_elliptic(f, q, a):
    """
    On input of a polynomial f that is a weil polynomial of degree 2 and has constant
    term q^a we check if it actually comes from an elliptic curve over GF(q^a).
    This uses theorem 4.1 of http://archive.numdam.org/article/ASENS_1969_4_2_4_521_0.pdf
    """
    if f[1] % q != 0:
        return True

    if a % 2 == 0:
        if (
            f[1] in [-2 * q ** (a // 2), 2 * q ** (a // 2)]
            or (q % 3 != 1 and f[1] in [-(q ** (a // 2)), q ** (a // 2)])
            or (q % 4 != 1 and f[1] == 0)
        ):
            return True
    else:
        if q in [2, 3]:
            if f[1] in [-(q ** ((a + 1) // 2)), q ** ((a + 1) // 2)]:
                return True
        if f[1] == 0:
            return True

    return False


def get_weil_polys(F):
    """
    Returns all degree 2 weil polynomials over F that are actually coming from an elliptic curve.
    """
    q = F.characteristic()
    a = F.degree()
    weil_polys = R.weil_polynomials(2, q**a)
    return [f for f in weil_polys if weil_polynomial_is_elliptic(f, q, a)]


def get_ordinary_weil_polys_from_values(q, a):
    weil_polys = R.weil_polynomials(2, q**a)
    return [f for f in weil_polys if f[1] % q != 0]


def eps_exp(alpha, eps, Sigma):
    return prod([sigma(alpha) ** my_pow for my_pow, sigma in zip(eps, Sigma)])


def gal_act_eps(eps, sigma):
    return tuple(eps[i - 1] for i in sigma)


def galois_action_on_embeddings(G_K):
    K = G_K.number_field()
    Kgal = G_K.splitting_field()
    embeddings = K.embeddings(Kgal)
    # first a shortcut in the case where G_K is normal in S_d since then it doesn't
    # matter for our application since we only care about the image
    # of the galois group in S_d
    d = K.absolute_degree()
    G_K_roots = TransitiveGroup(d, G_K.transitive_number())
    if G_K_roots.is_normal(SymmetricGroup(d)):
        id_G_K = G_K.Hom(G_K).identity()
        return G_K, id_G_K, id_G_K, Kgal, embeddings

    # now for the actual computation
    permutations = []
    for g in G_K.gens():
        phi = g.as_hom()
        g_perm = Permutation(
            [embeddings.index(phi * emb) + 1 for emb in embeddings]
        ).inverse()
        permutations.append(g_perm)
    G_K_emb = PermutationGroup(permutations, canonicalize=False)
    to_emb = G_K.hom(G_K_emb.gens())
    from_emb = G_K_emb.hom(G_K.gens())
    return G_K_emb, to_emb, from_emb, Kgal, embeddings


def class_group_as_additive_abelian_group(C_K):
    A = AdditiveAbelianGroup(C_K.gens_orders())

    def to_A(g):
        assert g in C_K
        return A(g.exponents())

    def from_A(a):
        assert a in A
        return C_K.prod(gi**ei for gi, ei in zip(C_K.gens(), a))

    return A, to_A, from_A


def ideal_push_forward(f, I):
    K = f.codomain()
    return K.ideal([f(i) for i in I.gens()])


def _class_group_norm_map_internal(phi):
    """
    If phi: L -> K is an embedding of number fields then this function returns
    the norm map from Cl_K to Cl_L induced by the relative norm from K to L.
    """
    L = phi.domain()
    K = phi.codomain()
    C_L = L.class_group()
    # currently only needed cause sage is really bad in doing
    # group morphisms of class groups
    A_L, to_A_L, from_A_L = class_group_as_additive_abelian_group(C_L)
    C_K = K.class_group()
    A_K, to_A_K, from_A_K = class_group_as_additive_abelian_group(C_K)

    Krel = K.relativize(phi, "z")
    to_K, from_K = Krel.structure()
    images = []
    for I in C_K.gens_ideals():
        I = ideal_push_forward(from_K, I)
        Nm_I = I.relative_norm()
        images.append(to_A_L(C_L(Nm_I)))
    norm_hom = A_K.hom(images)

    return norm_hom, to_A_K, from_A_L


def class_group_norm_map(phi, to_C_L=True):
    """
    If phi: L -> K is an embedding of number fields then this function returns
    the norm map from Cl_K to Cl_L induced by the relative norm from K to L.
    """
    norm_hom, to_A_K, from_A_L = _class_group_norm_map_internal(phi)
    C_K = phi.codomain().class_group()
    if to_C_L:

        def norm_map(I):
            I = C_K(I)
            return from_A_L(norm_hom(to_A_K(I)))

    else:

        def norm_map(I):
            I = C_K(I)
            return norm_hom(to_A_K(I))

    return norm_map


def split_embeddings(phi, embeddings):
    L = phi.domain()
    assert L.absolute_degree() == 2
    phi0 = phi(L.gen(0))
    im0 = embeddings[0](phi0)
    split1 = []
    split2 = []
    for embedding in embeddings:
        if embedding(phi0) == im0:
            split1.append(embedding)
        else:
            split2.append(embedding)
    return set(split1), set(split2)


def split_primes_iter(K, bound=oo, cache=True):
    if cache and not hasattr(K, "_prime_factorizations_cache"):
        K._prime_factorizations_cache = {}

    for p in primes(1, bound):
        if cache and p in K._prime_factorizations_cache:
            F = K._prime_factorizations_cache[p]
        else:
            F = (K * p).factor()
            if cache:
                K._prime_factorizations_cache[p] = F

        if not len(F) == K.absolute_degree():
            continue

        for pp, _ in F:
            yield pp


def primes_iter(K, bound=oo, cache=True):
    if cache and not hasattr(K, "_prime_factorizations_cache"):
        K._prime_factorizations_cache = {}

    for p in primes(1, bound):
        if cache and p in K._prime_factorizations_cache:
            F = K._prime_factorizations_cache[p]
        else:
            F = (K * p).factor()
            if cache:
                K._prime_factorizations_cache[p] = F

        for pp, _ in F:
            if pp.absolute_norm() < bound:
                yield pp


def one_aux_gen_list(C_K, class_group_gens, it):
    """Compute one Gen set"""
    running_class_group_gens = class_group_gens.copy()
    gen_list = []
    while running_class_group_gens:
        candidate = next(it)
        if candidate.smallest_integer() == 2:
            continue
        candidate_class = C_K(candidate)
        if candidate_class in running_class_group_gens:
            gen_list.append(candidate)
            running_class_group_gens.remove(candidate_class)
    return gen_list


def auxgens(K, auxgen_count=5):
    """Compute a list AuxGen of Gen sets"""

    C_K = K.class_group()
    class_group_gens = list(C_K.gens())

    it = primes_iter(K)

    aux_gen_list = [
        one_aux_gen_list(C_K, class_group_gens, it) for _ in range(auxgen_count)
    ]

    return aux_gen_list


def get_eps_type(eps):
    """Returns the type of an epsilon (quadratic, quartic, sextic), where
    an epsilon is considered as a tuple
    """

    if 6 in eps:
        if any(t in eps for t in [4, 8]):
            return "mixed"
        if len(set(eps)) == 1:
            # means it's all 6s
            return "type-2"
        return "quartic-non-constant"
    if any(t in eps for t in [4, 8]):
        if len(set(eps)) == 1:
            # means it's all 4s or all 8s
            return "sextic-constant"
        return "sextic-non-constant"
    if len(set(eps)) == 1:
        return "type-1"
    return "quadratic-non-constant"


def filter_ABC_primes(K, prime_list, eps_type):
    """Apply congruence and splitting conditions to primes in prime
    list, depending on the type of epsilon

    Args:
        K ([NumberField]): our number field
        prime_list ([list]): list of primes to filter
        eps_type ([str]): one of the possible epsilon types
    """

    if eps_type == "type-1":
        # no conditions
        return prime_list

    if eps_type == "quadratic-non-constant":
        # prime must split or ramify in K
        output_list = []

        for p in prime_list:
            if not K.ideal(p).is_prime():
                output_list.append(p)
        return output_list

    if eps_type == "sextic-non-constant":
        # prime must split or ramify in K, and be congruent to 2 mod 3
        output_list = []

        for p in prime_list:
            if p % 3 == 2:
                if not K.ideal(p).is_prime():
                    output_list.append(p)
        return output_list

    if eps_type == "sextic-constant":
        # prime must be congruent to 2 mod 3
        output_list = []

        for p in prime_list:
            if p % 3 == 2:
                output_list.append(p)
        return output_list

    if eps_type == "type-2":
        # prime must be congruent to 3 mod 4
        output_list = []

        for p in prime_list:
            if p % 4 == 3:
                output_list.append(p)
        return output_list

    if eps_type == "quartic-non-constant":
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
            if p % 12 == 11:
                if not K.ideal(p).is_prime():
                    output_list.append(p)
        return output_list

    raise ValueError("given type {} not a vaid epsilon type".format(eps_type))
