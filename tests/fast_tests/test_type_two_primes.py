import pytest
from sage_code.type_two_primes import get_type_2_primes, get_type_2_not_momose
from sage.all import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen

TEST_SETTINGS = {
    "bound": 1000,
}

x = polygen(QQ)
test_cases = [
    (x ** 2 + 1, 1),
    (x ** 3 + 3 * x ** 2 + 1, 1),
    (x ** 4 - 9 * x ** 2 - 10 * x + 50, 5),
]


@pytest.mark.parametrize("f, class_number", test_cases)
def test_get_type_2_primes(f, class_number):
    K = NumberField(f, "a")
    Kgal = K.galois_closure("b")
    embeddings = K.embeddings(Kgal)

    type_two_not_momose = get_type_2_not_momose(K, embeddings)
    _ = get_type_2_primes(K, embeddings, **TEST_SETTINGS)

    if class_number == 1:
        assert type_two_not_momose == set()
