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
