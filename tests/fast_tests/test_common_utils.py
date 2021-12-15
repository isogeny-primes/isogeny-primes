from sage.all import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring import polygen
from sage_code.common_utils import galois_action_on_embeddings


def test_galois_action_on_embeddings():
    x = polygen(QQ)
    # https://www.lmfdb.org/NumberField/5.1.35152.1
    f = x ** 5 - x ** 4 + 2 * x ** 3 - 4 * x ** 2 + x - 1
    K = NumberField(f, "a")
    G_K = K.galois_group()
    G_K_emb, to_emb, from_emb, Kgal, embeddings = galois_action_on_embeddings(G_K)
    # test that to_emb and from_emb are isomorphism
    assert G_K.hom(G_K) == from_emb * to_emb
    assert G_K_emb.hom(G_K_emb) == to_emb * from_emb
    assert to_emb.domain() == G_K
    assert to_emb.codomain() == G_K_emb
    assert from_emb.domain() == G_K_emb
    assert from_emb.codomain() == G_K
    assert embeddings[0] == G_K.gen(0).as_hom() * embeddings[G_K_emb.gen(0)(1) - 1]
    assert Kgal.degree() == G_K_emb.cardinality() == 20
    assert G_K_emb.degree() == len(embeddings) == 5
