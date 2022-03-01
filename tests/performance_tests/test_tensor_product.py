"""test_tensor_product.py

Tests for tensor products.

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

from sage.all import QQ
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.special import block_diagonal_matrix

M = MatrixSpace(QQ, 2)
one = M.one()
matrices = [M.random_element() for i in range(2000)]


def test_m_tensor_one():
    for m in matrices:
        m.tensor_product(one, subdivide=False)


def test_one_tensor_m():
    for m in matrices:
        one.tensor_product(m, subdivide=False)


def test_m_tensor_one_subdivide():
    for m in matrices:
        m.tensor_product(one, subdivide=True)


def test_one_tensor_m_subdivide():
    for m in matrices:
        one.tensor_product(m, subdivide=True)


def test_one_tensor_m_block_diag():
    for m in matrices:
        block_diagonal_matrix([m] * one.nrows())
