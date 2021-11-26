"""common_utils.py

    Common methods and classes for multiple parts of the routine.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    ====================================================================

"""

from sage.all import PolynomialRing, Rationals

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
    weil_polys = R.weil_polynomials(2, q ** a)
    return [f for f in weil_polys if weil_polynomial_is_elliptic(f, q, a)]
