"""test_common_utils.py

Tests for the common utils functions.

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
from sage.all import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen
from sage_code.common_utils import galois_action_on_embeddings, is_b_smooth

x = polygen(QQ)
test_cases = [
    (x**5 - x**4 + 2 * x**3 - 4 * x**2 + x - 1, 20),
    (x**6 - 2 * x**5 + 6 * x**4 + 22 * x**3 + 41 * x**2 + 48 * x + 36, 6),
]


@pytest.mark.parametrize("f, gal_deg", test_cases)
def test_galois_action_on_embeddings(f, gal_deg):
    K = NumberField(f, "a")
    G_K = K.galois_group()
    G_K_emb, to_emb, from_emb, Kgal, embeddings = galois_action_on_embeddings(G_K)
    # test that to_emb and from_emb are isomorphism
    assert G_K.hom(G_K) == from_emb * to_emb
    assert G_K_emb.hom(G_K_emb) == to_emb * from_emb
    assert to_emb.domain() == G_K
    assert to_emb.codomain() == G_K_emb
    assert from_emb.domain() == G_K_emb
    assert from_emb.codomain() == G_K
    assert embeddings[0] == G_K.gen(0).as_hom() * embeddings[G_K_emb.gen(0)(1) - 1]
    assert Kgal.degree() == G_K_emb.cardinality() == gal_deg
    assert G_K_emb.degree() == len(embeddings) == K.degree()


@pytest.mark.parametrize(
    "n, b, expected",
    [
        [1, 3, (True, [])],
        [2 * 3 * 37 * 41, 17, (False, [2, 3, 37 * 41])],
        [2 * 3 * 37 * 41, 40, (False, [2, 3, 37, 41])],
        [2 * 3 * 37 * 41, 42, (True, [2, 3, 37, 41])],
    ],
)
def test_is_b_smooth(n, b, expected):
    result = is_b_smooth(n, b)
    assert result == expected
