"""test_is_coprime.py

Tests for coprimality.

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

from sage.all import QQ, prime_range
from sage.arith.misc import next_prime, gcd
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen

x = polygen(QQ)
K = NumberField(x ** 2 - 11 * 17 * 9011 * 23629, "b")
test_cases_true = [(K, int(q), q, False) for q in prime_range(1000)]
test_cases_false = [(K, int(q), next_prime(q), True) for q in prime_range(1000)]


def test_warmup():
    for K, q, p, result in test_cases_true + test_cases_false:
        qq = (q * K).factor()[0][0]
        assert qq.is_coprime(p) == (not qq.divides(p)) == result


def test_is_coprime():
    for K, q, p, result in test_cases_true + test_cases_false:
        qq = (q * K).factor()[0][0]
        is_coprime = qq.is_coprime(p)
        assert is_coprime == result


def test_divides():
    for K, q, p, result in test_cases_true + test_cases_false:
        qq = (q * K).factor()[0][0]
        is_coprime = not qq.divides(p)
        assert is_coprime == result


def test_norm_gcd():
    for K, q, p, result in test_cases_true + test_cases_false:
        qq = (q * K).factor()[0][0]
        is_coprime = gcd(qq.absolute_norm(), p) == 1
        assert is_coprime == result


def test_norm_div():
    for K, q, p, result in test_cases_true + test_cases_false:
        qq = (q * K).factor()[0][0]
        is_coprime = qq.absolute_norm() % p != 0
        assert is_coprime == result


def test_norm_divides():
    for K, q, p, result in test_cases_true + test_cases_false:
        qq = (q * K).factor()[0][0]
        is_coprime = not p.divides(qq.absolute_norm())
        assert is_coprime == result


def test_smallest_int_eq():
    for K, q, p, result in test_cases_true + test_cases_false:
        qq = (q * K).factor()[0][0]
        is_coprime = qq.smallest_integer() != p
        assert is_coprime == result
