from sage.misc.misc_c import prod
from sage.rings.number_field.number_field import QuadraticField

from sage_code import ideal_log_relation


def test_ideal_log_relation():
    K = QuadraticField(-41)
    q = (K * 3).factor()[0][0]
    C_K = K.class_group()
    t = ideal_log_relation(q)
    e = q.ideal_class_log()
    assert q == t * prod(gi ** ei.sage() for gi, ei in zip(C_K.gens_ideals(), e))
