"""
This init file contains code that monkey patches sage so that :trac:`32910` is fixed::
The entire content of this file can be made empty once this is fixed in upstream sage
"""
import sage
from sage.all import ZZ, diagonal_matrix, vector, matrix


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
