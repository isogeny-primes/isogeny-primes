import pytest

from sage.all import (
    QuadraticField,
    companion_matrix,
    prime_range,
)  # pylint: disable=no-name-in-module
from sage_code.weeding import (
    get_dirichlet_character,
    works_method_of_appendix,
    is_rank_of_twist_zero,
)


def test_get_dirichlet_character():
    K = QuadraticField(-31)
    chi = get_dirichlet_character(K)
    assert (3 * K).is_prime()
    assert chi(3) == -1
    assert not (5 * K).is_prime()
    assert chi(5) == 1
    assert (73 * K).is_prime()
    assert chi(73) == -1


@pytest.mark.parametrize(
    "D, p, works",
    [
        (5, 79, True),
        (-31, 73, False),
    ],
)
def test_works_method_of_appendix(D, p, works):
    K = QuadraticField(D)
    result = works_method_of_appendix(p, K)
    assert result == works


def test_is_rank_of_twist_zero():
    p = 73
    K = QuadraticField(-31)
    chi = get_dirichlet_character(K)
    assert not is_rank_of_twist_zero(p, chi)
