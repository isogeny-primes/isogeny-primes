"""frobenius_polynomials.py

Code used to generate isogenies of given signatures.

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

from sage.all import ZZ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.polynomial.polynomial_ring import polygen
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_number_field import EllipticCurve_number_field
from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
from sage.rings.integer import Integer


def uniformizer(p: NumberFieldFractionalIdeal):
    for g in p.gens_two():
        if g.valuation(p) == 1:
            return g
    raise ValueError(f"No uniformizer of {p} found.")


# This is needed cause the reduction function of sage has errors
def reduction(E: EllipticCurve_number_field, p: NumberFieldFractionalIdeal):
    if not E.has_good_reduction(p):
        raise ValueError(f"The curve must have good reduction at the place.")
    Ep = E.local_data(p).minimal_model()
    K = E.base_ring()
    Fp = K.residue_field(p)
    Fp2 = GF(Fp.order(), "abar", modulus=Fp.modulus())
    iso = Fp.hom([Fp2.gen(0)], Fp2)
    Ebar = EllipticCurve([iso(Fp(ai)) for ai in Ep.ainvs()])
    return Ebar


kodaira_to_e_dict = {
    "II": 6,
    "III": 4,
    "IV": 3,
    "II*": 6,
    "III*": 4,
    "IV*": 3,
}


def semistable_ramification(local_data, p=None):
    if p:
        local_data = local_data.local_data(p)
    if local_data.has_good_reduction():
        return 1
    if local_data.has_multiplicative_reduction():
        return 1
    kodaira = str(local_data.kodaira_symbol())
    try:
        return kodaira_to_e_dict[kodaira]
    except KeyError:
        return 2


def semi_stable_frobenius_polynomial(
    E: EllipticCurve_number_field, q: NumberFieldFractionalIdeal, t=1
):
    """
    Input:
        E - an elliptic curve over a number field K
        q - a prime of the number field K
        t - an element of t coprime to q, used to change the purely
            ramified extension that is chosen

    Output:
        A frobenius polymial of an elliptic curve E' over O_K/q that
        is obtained by computing a semistable model of E over a purely
        ramified extension of K_q.
        Note that in the case of multiplicative reduction we return:
        (x-1)*(x-Norm(q)) in the split case and
        (x+1)*(x+Norm(q)) in the nonsplit case

    Note that the output is slightly random, since it is dependent ony the
    purely ramified extension that is chosen. So E' is and hence the frobenius
    polynomial is only defined up to twists.
    """
    K = E.base_ring()
    assert not K(t) in q

    local_data = E.local_data(q)

    if local_data.has_good_reduction():
        Ebar = E.reduction(q)
    elif E.j_invariant().valuation(q) >= 0:
        x = polygen(K)
        u = uniformizer(q)
        e = semistable_ramification(local_data)
        L = K.extension(x ** e - t * u, "a")
        EL = E.change_ring(L)
        qL = L.primes_above(q)[0]
        Ebar = reduction(EL, qL)
    else:
        # we have potentially multiplicative reduction
        x = polygen(ZZ)
        if local_data.has_nonsplit_multiplicative_reduction():
            return (x + 1) * (x + q.norm())
        else:
            return (x - 1) * (x - q.norm())

    return Ebar.frobenius_polynomial()


def isogeny_character_values(
    E: EllipticCurve_number_field, p: Integer, q: NumberFieldFractionalIdeal
):
    f = semi_stable_frobenius_polynomial(E, q)
    roots = f.roots(GF(p), multiplicities=False)
    if not roots:
        raise ValueError(
            f"The (semi stable model of) the curve does not admit an {p}-isogeny over F_q"
        )
    else:
        return roots


def isogeny_character_values_12(
    E: EllipticCurve_number_field, p: Integer, q: NumberFieldFractionalIdeal
):
    values = isogeny_character_values(E, p, q)
    return [a ** 12 for a in values]
