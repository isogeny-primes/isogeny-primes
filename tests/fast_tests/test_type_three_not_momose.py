import pytest
from sage_code.type_three_not_momose import type_three_not_momose
from sage_code.pre_type_one_two import get_strong_type_3_epsilons
from sage.all import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen

TEST_SETTINGS = {
    "bound": 1000,
}

x = polygen(QQ)
test_cases = [
    (x ** 2 + 1, 1, 1),
    (x ** 3 + 3 * x ** 2 + 1, 1, 0),
    (x ** 4 - 9 * x ** 2 - 10 * x + 50, 5, 1),
    (x ** 4 - x ** 3 - x ** 2 - 2 * x + 4, 1, 2),
]


@pytest.mark.parametrize("f, class_number, strong_L_count", test_cases)
def test_type_three_not_momose(f, class_number, strong_L_count):
    K = NumberField(f, "a")
    Kgal = K.galois_closure("b")
    embeddings = K.embeddings(Kgal)

    strong_type_3_epsilons = get_strong_type_3_epsilons(K, embeddings)

    assert len(strong_type_3_epsilons) == 2 * strong_L_count

    type_3_not_momose_primes, _ = type_three_not_momose(
        K, embeddings, strong_type_3_epsilons
    )

    if class_number == 1:
        assert type_3_not_momose_primes == []
