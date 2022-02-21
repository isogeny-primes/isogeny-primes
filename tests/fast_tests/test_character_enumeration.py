"""test_character_enumeration.py

Tests for the ICE filter.

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

import pytest
from sage.all import GF, NumberField, QQ
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring import polygen

from sage_code.character_enumeration import (
    filter_possible_values,
    lifts_in_hasse_range,
    get_prime_gens,
)


@pytest.mark.parametrize(
    "possible_values_list, q, e, prime, result",
    [
        # +/-1 and +/-q^e should always be allowed
        ([1, 2, -1, -2], 2, 1, 271, [1, 2, -1, -2]),
        ([1, 4, -1, -4], 2, 2, 271, [1, 4, -1, -4]),
        # test_isogeny_character_values_12 shows these are valid
        ([13, 57], 11, 1, 73, [13, 57]),
    ],
)
def test_filter_possible_values(possible_values_list, q, e, prime, result):
    prime_field = GF(prime)
    possible_values_list = [prime_field(v) for v in possible_values_list]
    possible_values = filter_possible_values(possible_values_list, q, e, prime_field)
    assert possible_values == result


@pytest.mark.parametrize(
    "fq, res_class, expected_range",
    [
        (11, GF(3)(0), [0, -3, -6, 3, 6]),
        # make sure we include +/-2sqrt(fq) but not +/-(2sqrt(fq)+1),
        (145757 ** 2, GF(357686312646216567629137)(2 * 145757), [2 * 145757]),
        (145757 ** 2, GF(357686312646216567629137)(-2 * 145757), [-2 * 145757]),
        (145757 ** 2, GF(357686312646216567629137)(2 * 145757 + 1), []),
        (145757 ** 2, GF(357686312646216567629137)(-2 * 145757 - 1), []),
    ],
)
def test_lifts_in_hasse_range(fq, res_class, expected_range):
    result = lifts_in_hasse_range(fq, res_class)
    assert result == expected_range


x = polygen(QQ)


@pytest.mark.parametrize(
    "f, q",
    [
        (x ** 2 - x + 1007, 13),
        (x ** 2 - x + 1007, 17),
        (x ** 2 - x + 1007, 19),
        (x ** 2 - 11 * 17 * 9011 * 23629, 17),
        (x ** 2 - 11 * 17 * 9011 * 23629, 47),
        (x ** 2 - 11 * 17 * 9011 * 23629, 89),
    ],
)
@pytest.mark.skip(reason="Only for profiling putative faster solution")
def test_ideal_log_relation_prime_gens(f, q):
    # also tests get_prime_gens
    # x = polygen(QQ)
    # f = x ** 2 - x + 1007
    # q = 13
    K = NumberField(f, "b")
    C_K = K.class_group()
    qq = (q * K).factor()[0][0]

    prime_gens, relations = get_prime_gens(C_K)
    ideal_gens = C_K.gens_ideals()
    for prime, ideal, relation in zip(prime_gens, ideal_gens, relations):
        assert prime * relation == ideal

    exponents = C_K(qq).exponents()

    t = qq.ideal_log_relation()
    alpha_new = t * prod([(relation ** e) for relation, e in zip(relations, exponents)])

    sanity_check = prod([Q ** a for Q, a in zip(prime_gens, exponents)])
    assert C_K(sanity_check) == C_K(qq)

    the_principal_ideal = qq * prod([Q ** (-a) for Q, a in zip(prime_gens, exponents)])

    alphas = the_principal_ideal.gens_reduced()
    assert len(alphas) == 1, "the principal ideal isn't principal!!!"
    alpha = alphas[0]

    assert alpha == alpha_new or (alpha / alpha_new).is_unit()
