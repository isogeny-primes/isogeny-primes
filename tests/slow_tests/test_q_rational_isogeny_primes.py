from sage.all import QQ, polygen

from isogeny_primes import get_isogeny_primes
from sage_code.common_utils import EC_Q_ISOGENY_PRIMES
import logging


logger = logging.getLogger(__name__)

TEST_SETTINGS = {
    "norm_bound": 20,
    "bound": 1000,
    "loop_curves": False,
    "heavy_filter": True,
}


def test_rational_isogeny_primes():
    x = polygen(QQ)
    K = QQ.extension(x - 1, "D")

    superset = get_isogeny_primes(K, **TEST_SETTINGS)
    improperly_ruled_out = EC_Q_ISOGENY_PRIMES.difference(superset)
    assert improperly_ruled_out == set()
    todo = set(superset).difference(EC_Q_ISOGENY_PRIMES)
    # would be nice if we could automatically do QQ
    # i.e. we could rule out 23 as well.
    assert todo == set([23])
