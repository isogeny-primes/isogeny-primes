"""
Make sure we show a decent error if a to old version of sage is being used
"""
from sage.structure.proof.proof import get_flag

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


def _cache_bnfisprincipal(self, proof=None, gens=False):
    r"""
    This function is essentially the implementation of
    :meth:`is_principal`, :meth:`gens_reduced` and
    :meth:`ideal_class_log`.
    INPUT:
    - ``self`` -- an ideal
    - ``proof`` -- proof flag.  If ``proof=False``, assume GRH.
    - ``gens`` -- (default: False) if True, also computes the reduced
      generators of the ideal.
    OUTPUT:
    None.  This function simply caches the results: it sets
    ``_ideal_class_log`` (see :meth:`ideal_class_log`),
    ``_is_principal`` (see :meth:`is_principal`) and
    ``_reduced_generators``.
    TESTS:
    Check that no warnings are triggered from PARI/GP (see :trac:`30801`)::
        sage: K.<a> = NumberField(x^2 - x + 112941801)
        sage: I = K.ideal((112941823, a + 49942513))
        sage: I.is_principal()
        False
    """
    # Since pari_bnf() is cached, this call to pari_bnf() should not
    # influence the run-time much.  Also, this simplifies the handling
    # of the proof flag: if we computed bnfisprincipal() in the past
    # with proof=False, then we do not need to recompute the result.
    # We just need to check correctness of pari_bnf().
    proof = get_flag(proof, "number_field")
    bnf = self.number_field().pari_bnf(proof)

    # If we already have _reduced_generators, no need to compute them again
    if hasattr(self, "_reduced_generators"):
        gens = False

    # Is there something to do?
    if hasattr(self, "_ideal_class_log") and not gens:
        self._is_principal = not any(self._ideal_class_log)
        return

    if not gens:
        v = bnf.bnfisprincipal(self.pari_hnf(), 0)
        self._ideal_class_log = list(v)
        self._is_principal = not any(self._ideal_class_log)
    else:
        # TODO: this is a bit of a waste. We ask bnfisprincipal to compute the compact form and then
        # convert this compact form back into an expanded form.
        # (though calling with 3 instead of 5 most likely triggers an error with memory allocation failure)
        v = bnf.bnfisprincipal(self.pari_hnf(), 5)
        e = v[0]
        t = v[1]
        t = bnf.nfbasistoalg(bnf.nffactorback(t))
        self._ideal_class_log = list(e)
        g = self.number_field()(t)
        self._ideal_log_relation = g
        self._is_principal = not any(self._ideal_class_log)

        if self._is_principal:
            self._reduced_generators = (g,)
        elif gens:
            # Non-principal ideal
            self._reduced_generators = self.gens_two()


def ideal_log_relation(self, proof=None):
    if not hasattr(self, "_ideal_log_relation"):
        self._cache_bnfisprincipal(proof=proof, gens=True)

    return self._ideal_log_relation


sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal._cache_bnfisprincipal = (
    _cache_bnfisprincipal
)
sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal.ideal_log_relation = (
    ideal_log_relation
)
