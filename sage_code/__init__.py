"""
Make sure we show a decent error if a to old version of sage is being used
"""
try:
    from sagemath.check_version import check_version

    check_version(">=9.4")
except ModuleNotFoundError as err:
    raise RuntimeError(
        """
        Not all requirements of Isogeny Primes are installed. Please do
            sage -pip install -r requirements.txt
        Before continuing.
        """
    ) from ModuleNotFoundError

"""
The code in the rest of this file monkey patches sage so that :trac:`32910` is fixed::
This can be removed once this is fixed in upstream sage
"""

import sage
from sage.all import ZZ, IntegerModRing, gcd, prod, diagonal_matrix, vector, matrix
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement


def __contains__(self, x):
    """
    Test whether ``x`` is an element of this subgroup.

    EXAMPLES::

        sage: G.<a,b> = AbelianGroup(2)
        sage: A = G.subgroup([a])
        sage: a in G
        True
        sage: a in A
        True

    TESTS:

    Check that :trac:`32910` is fixed::

        sage: G.<a,b> = AbelianGroup(2, [4, 576])
        sage: Hgens = [a^2, a*b^2]
        sage: H = G.subgroup(Hgens)
        sage: [g in H for g in (a^3, b^2, b^3, a^3*b^2, "junk")]
        [False, False, False, True, False]

    Check that :trac:`31507` is fixed::

        sage: G = AbelianGroup(2, gens_orders=[16, 16])
        sage: f0, f1 = G.gens()
        sage: H = G.subgroup([f0*f1^3])
        sage: [g in H for g in (f0, f0*f1^2, f0*f1^3, f0*f1^4)]
        [False, False, True, False]

        sage: G.<a,b> = AbelianGroup(2)
        sage: Hgens =  [a*b, a*b^-1]
        sage: H = G.subgroup(Hgens)
        sage: b^2 in H
        True
    """

    if not isinstance(x, AbelianGroupElement):
        return False
    if x.parent() is self:
        return True
    elif x in self.ambient_group():
        amb_inv = self.ambient_group().gens_orders()
        inv_basis = diagonal_matrix(ZZ, amb_inv)
        gens_basis = matrix(
            ZZ, len(self._gens), len(amb_inv), [g.list() for g in self._gens]
        )
        return vector(ZZ, x.list()) in inv_basis.stack(gens_basis).row_module()
    return False


sage.groups.abelian_gps.abelian_group.AbelianGroup_subgroup.__contains__ = __contains__

# the following is to significantly speed up the calculation of the dirichlet
# group of a number field.
def _splitting_classes_gens_(K, m, d):
    r"""
    Given a number field `K` of conductor `m` and degree `d`,
    this returns a set of multiplicative generators of the
    subgroup of `(\mathbb{Z}/m\mathbb{Z})^{\times}/(\mathbb{Z}/m\mathbb{Z})^{\times d}`
    containing exactly the classes that contain the primes splitting
    completely in `K`.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field import _splitting_classes_gens_
        sage: K = CyclotomicField(101)
        sage: L = K.subfields(20)[0][0]
        sage: L.conductor()
        101
        sage: _splitting_classes_gens_(L,101,20)
        [95]

        sage: K = CyclotomicField(44)
        sage: L = K.subfields(4)[0][0]
        sage: _splitting_classes_gens_(L,44,4)
        [37]

        sage: K = CyclotomicField(44)
        sage: L = K.subfields(5)[0][0]
        sage: K.degree()
        20
        sage: L
        Number Field in zeta44_0 with defining polynomial x^5 - 2*x^4 - 16*x^3 + 24*x^2 + 48*x - 32 with zeta44_0 = 3.837971894457990?
        sage: L.conductor()
        11
        sage: _splitting_classes_gens_(L,11,5)
        [10]

    """

    R = K.ring_of_integers()
    Zm = IntegerModRing(m)
    unit_gens = Zm.unit_gens()
    ZZunits = ZZ ** len(unit_gens)
    unit_relations = [gcd(d, x.multiplicative_order()) for x in unit_gens]
    # sparse=False can be removed if the code below doesn't raise the following
    # AttributeError: 'Matrix_integer_sparse' object has no attribute '_clear_denom'
    D = diagonal_matrix(unit_relations, sparse=False)
    Zmstar = ZZunits / D.row_module()

    def map_Zmstar_to_Zm(h):
        li = h.lift().list()
        return prod(unit_gens[i] ** li[i] for i in range(len(unit_gens)))

    Hgens = []
    H = Zmstar.submodule([])
    Horder = Zmstar.cardinality() / d

    for g in Zmstar:
        if H.cardinality() == Horder:
            break
        if g not in H:
            u = map_Zmstar_to_Zm(g)
            p = u.lift()
            while not p.is_prime():
                p += m
            f = R.ideal(p).prime_factors()[0].residue_class_degree()
            h = g * f
            if h not in H:
                Hgens += [h]
                H = Zmstar.submodule(Hgens)

    return [map_Zmstar_to_Zm(h) for h in Hgens]


sage.rings.number_field.number_field._splitting_classes_gens_ = _splitting_classes_gens_
