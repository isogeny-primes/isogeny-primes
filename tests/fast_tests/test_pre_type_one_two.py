import pytest
from sage.arith.misc import is_squarefree
from sage.rings.number_field.number_field import QuadraticField

from sage_code.pre_type_one_two import contains_imaginary_quadratic_field

square_free_D = [D for D in range(-100, 100) if is_squarefree(D) and D != 1]


@pytest.mark.parametrize("D", square_free_D)
def test_contains_imaginary_quadratic_field(D):
    K = QuadraticField(D)
    result = contains_imaginary_quadratic_field(K)
    assert result == (D < 0)
