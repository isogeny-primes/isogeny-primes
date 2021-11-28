import pytest
from sage.all import GF

from sage_code.character_enumeration import filter_possible_values


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
