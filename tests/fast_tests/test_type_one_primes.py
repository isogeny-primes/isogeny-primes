import pytest

from sage_code.type_one_primes import (
    is_formall_immersion,
    cached_bad_formal_immersion_data,
)


@pytest.mark.parametrize(
    "d,p,expected_result",
    [
        [3, 71, 1],
        [4, 71, 1],
        [5, 71, 1],
        [6, 71, 1],
        [7, 71, 0],
    ],
)
def test_is_formall_immersion(d, p, expected_result):
    assert is_formall_immersion(d, p) == expected_result


@pytest.mark.parametrize("d", range(1, 13))
def test_cached_bad_formal_immersion_data(d):
    cached_bad_formal_immersion_data(d)
