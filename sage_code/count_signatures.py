"""count_signatures.py

    Combinatorial counts of isogeny signatures, used in the section of the
    paper that quantifies how many signatures remain after the main theorem.

    A *signature* of degree ``d`` is an element of ``{0,4,6,8,12}^d``; its
    *trace* is the sum of its entries.  Writing ``K`` for the Galois closure of
    a degree ``d`` field ``k``, the group ``G x H`` (with ``G = Gal(K/Q)`` and
    ``H = Aut(k)``) acts on signatures, and we also identify a signature with
    its dual ``12 - eps``.  The quantities computed here are:

      * ``N(d)``  -- the minimal number of signature classes in
        ``(G x H)\\Sigma_k`` over all degree ``d`` fields ``k``.  The minimum
        is attained exactly when ``Gal(K/Q) = S_d`` (the generic case), in
        which signature classes correspond to content vectors up to reversal.

      * the number of *remaining* classes, i.e. those the main theorem does not
        resolve: trace ``= 0 (mod 12)`` (the obstruction that remains even
        under GRH) and trace ``= 0 (mod 6)`` (the unconditional obstruction),
        excluding the single Type 1 class.

    This module is deliberately free of any Sage dependency so that it can be
    run directly to reproduce the table in the paper and to cross-check the
    closed-form formulas against brute-force enumeration:

        python3 sage_code/count_signatures.py

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2022 Barinder S. Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    The authors can be reached at: barinder.s.banwait@gmail.com and
    maarten@mderickx.nl.

    ====================================================================

"""

from itertools import product
from math import comb

# The possible values of a signature at a single embedding.
VALUES = (0, 4, 6, 8, 12)


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------
def N(d):
    """Minimal number of signature classes ``|(G x H)\\Sigma_k|`` over all
    degree ``d`` number fields ``k``.

    Equals the number of content vectors ``(n0,...,n4)`` with ``sum = d`` taken
    up to the reversal ``(n0,n1,n2,n3,n4) -> (n4,n3,n2,n1,n0)``.  Burnside gives

        N(d) = ( C(d+4, 4) + C(floor(d/2)+2, 2) ) / 2.
    """
    return (comb(d + 4, 4) + comb(d // 2 + 2, 2)) // 2


def raw_remaining_mod12(d):
    """Number of signatures (out of ``5^d``) with trace ``= 0 (mod 12)``.

    Closed form ``(5^d + 2*2^d + 3^d) / 6`` from the roots-of-unity filter.
    """
    return (5 ** d + 2 * 2 ** d + 3 ** d) // 6


def raw_remaining_mod6(d):
    """Number of signatures (out of ``5^d``) with trace ``= 0 (mod 6)``.

    Closed form ``(5^d + 2*2^d) / 3`` from the roots-of-unity filter.
    """
    return (5 ** d + 2 * 2 ** d) // 3


# ---------------------------------------------------------------------------
# Brute-force enumeration (used to check the closed forms)
# ---------------------------------------------------------------------------
def content_vectors(d):
    """Yield every ``(n0,n1,n2,n3,n4)`` of non-negative integers summing to d."""
    for n0 in range(d + 1):
        for n1 in range(d - n0 + 1):
            for n2 in range(d - n0 - n1 + 1):
                for n3 in range(d - n0 - n1 - n2 + 1):
                    yield (n0, n1, n2, n3, d - n0 - n1 - n2 - n3)


def _content_trace(content):
    n0, n1, n2, n3, n4 = content
    return 4 * n1 + 6 * n2 + 8 * n3 + 12 * n4


def generic_class_counts(d):
    """Counts in the generic case ``Gal(K/Q) = S_d``, where signature classes
    are content vectors up to reversal.

    Returns ``(total, remaining_mod12, remaining_mod6)`` where the remaining
    counts exclude the single Type 1 class (the all-0 and all-12 contents, which
    coincide under reversal, and which the main theorem resolves).
    """
    seen = set()
    total = rem12 = rem6 = 0
    type1 = {(d, 0, 0, 0, 0), (0, 0, 0, 0, d)}
    for content in content_vectors(d):
        if content in seen:
            continue
        seen.add(content)
        seen.add(content[::-1])
        total += 1
        if content in type1:
            continue
        trace = _content_trace(content)
        if trace % 12 == 0:
            rem12 += 1
        if trace % 6 == 0:
            rem6 += 1
    return total, rem12, rem6


def raw_class_counts(d):
    """Brute-force counts over all ``5^d`` signatures: ``(total, #trace=0 mod12,
    #trace=0 mod6)``."""
    total = rem12 = rem6 = 0
    for eps in product(VALUES, repeat=d):
        total += 1
        trace = sum(eps)
        if trace % 12 == 0:
            rem12 += 1
        if trace % 6 == 0:
            rem6 += 1
    return total, rem12, rem6


# ---------------------------------------------------------------------------
# Self-check and table
# ---------------------------------------------------------------------------
def _check(max_brute=10):
    """Assert the closed forms agree with brute force for d up to max_brute."""
    for d in range(1, max_brute + 1):
        total, _, _ = generic_class_counts(d)
        assert total == N(d), (d, total, N(d))
        _, raw12, raw6 = raw_class_counts(d)
        assert raw12 == raw_remaining_mod12(d), (d, raw12, raw_remaining_mod12(d))
        assert raw6 == raw_remaining_mod6(d), (d, raw6, raw_remaining_mod6(d))
    # Spot checks tying into the paper's existing numbers.
    assert N(2) == 9
    assert N(3) == 19
    assert raw_remaining_mod12(2) == 7  # 7/25 of the 25 quadratic signatures


def main():
    _check()
    print("Generic case  Gal(K/Q) = S_d  (signature classes = content up to reversal)")
    print(
        "{:>3}  {:>8}  {:>10}  {:>9}  {:>10}  {:>9}".format(
            "d", "N(d)", "rem(12)", "prop(12)", "rem(6)", "prop(6)"
        )
    )
    for d in range(2, 11):
        total, rem12, rem6 = generic_class_counts(d)
        print(
            "{:>3}  {:>8}  {:>10}  {:>9.4f}  {:>10}  {:>9.4f}".format(
                d, total, rem12, rem12 / total, rem6, rem6 / total
            )
        )
    print("\nlimits as d -> infinity:  prop(12) -> 1/6 = {:.4f}   prop(6) -> 1/3 = {:.4f}".format(1 / 6, 1 / 3))

    print("\nRaw count over all 5^d signatures (independent of the field k)")
    print(
        "{:>3}  {:>10}  {:>10}  {:>9}  {:>10}  {:>9}".format(
            "d", "5^d", "rem(12)", "prop(12)", "rem(6)", "prop(6)"
        )
    )
    for d in range(2, 11):
        tot = 5 ** d
        r12 = raw_remaining_mod12(d)
        r6 = raw_remaining_mod6(d)
        print(
            "{:>3}  {:>10}  {:>10}  {:>9.4f}  {:>10}  {:>9.4f}".format(
                d, tot, r12, r12 / tot, r6, r6 / tot
            )
        )


if __name__ == "__main__":
    main()
