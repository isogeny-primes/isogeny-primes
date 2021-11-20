import pytest
from sage_code.type_one_primes import is_formall_immersion


@pytest.mark.parametrize(
    "d,p,expected_result",
    [
        [3, 71, 1],
        [4, 71, 1],
        [5, 71, 1],
        [6, 71, 9],
        [7, 71, 0],
    ],
)
def test_is_formall_immersion(d, p, expected_result):
    assert is_formall_immersion(d, p) == expected_result
