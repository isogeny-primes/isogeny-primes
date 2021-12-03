import pytest
from sage.all import GF

from sage_code.character_enumeration import filter_possible_values, lifts_in_hasse_range


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
