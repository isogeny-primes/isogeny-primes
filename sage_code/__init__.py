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

        sage: Zmstar.<a,b> = AbelianGroup(2,[4,576])
        sage: Hgens =  [a**2,a*b**2]
        sage: H = Zmstar.subgroup(Hgens)
        sage: g = Zmstar.gen(1)**3
        sage: g in H
        False

    """
    amb_inv = self.ambient_group().gens_orders()
    inv_basis = diagonal_matrix(ZZ, amb_inv)
    gens_basis = matrix(
        ZZ, len(self._gens), len(amb_inv), [g.list() for g in self._gens]
    )
    return vector(ZZ, x.list()) in inv_basis.stack(gens_basis).row_module()


sage.groups.abelian_gps.abelian_group.AbelianGroup_subgroup.__contains__ = __contains__
