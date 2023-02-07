"""isogeny_primes.py

    Return finite list of isogeny primes attached to a number field.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2022 Barinder S. Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or any later version.

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

# Imports

import argparse
import logging

from sage.all import Integer, prime_divisors, is_prime

from sage_code.common_utils import is_b_smooth
from sage_code.strong_uniform_bounds import unif_bd, type_one_unif_bound


def do_uniform(d, trial_division_bound, aux_bound):

    epsilons = [
        tuple([trace] + [0] * (d - 1)) for trace in range(0, 6 * d, 2) if trace != 2 and (trace == 0 or trace % 6 != 0)
    ]

    for eps in epsilons:
        do_uniform_eps(d, eps, trial_division_bound, aux_bound, print_trace=True)


def do_uniform_type_1(d, trial_division_bound, aux_bound):
    mult_bound = Integer(type_one_unif_bound(d, aux_bound))
    is_smooth, factors = is_b_smooth(mult_bound, trial_division_bound)
    factors_str = "{" + ", ".join(str(i) for i in factors) + "}"
    logging.info(f"A superset of the type 1 isogeny primes (trace eps = 0 and {12*d}) is: {factors_str}")
    if not is_smooth and not is_prime(factors[-1]):
        logging.warning(f"unable to factor by trial division, the last factor is not prime!")


def do_uniform_eps(d, eps, trial_division_bound, aux_bound, print_trace=False):

    if eps == tuple(d * [0]):
        do_uniform_type_1(d, trial_division_bound, aux_bound)
        return

    mult_bound = Integer(unif_bd(d, eps, aux_bound))
    is_smooth, factors = is_b_smooth(mult_bound, trial_division_bound)
    factors_str = "{" + ", ".join(str(i) for i in factors) + "}"
    trace = sum(eps)
    if print_trace:
        info = f"A superset of the isogeny primes for trace eps = {trace} and {12*d - trace} is: {factors_str}"
    else:
        info = f"A superset of the isogeny primes for trace eps = trace {eps} = {trace} is: {factors_str}"
    logging.info(info)
    if not is_smooth and not is_prime(factors[-1]):
        logging.warning(f"unable to factor by trial division, the last factor is not prime!")


def cli_handler(args):  # pylint: disable=redefined-outer-name

    loglevel = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        level=loglevel,
    )
    logging.debug("Debugging level for log messages set.")

    if args.eps:
        if args.d < len(args.eps):
            raise ValueError(f"eps should have length at most d={args.d}, but has length {len(args.eps)} instead")
        eps = tuple(args.eps + [0] * (args.d - len(args.eps)))
        do_uniform_eps(args.d, eps, args.trial_division_bound, args.aux_bound)
        return

    do_uniform(args.d, args.trial_division_bound, args.aux_bound)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "d",
        metavar="d",
        type=int,
        help="the degree of the number field for which the bound is computed",
    )
    parser.add_argument(
        "--eps",
        required=False,
        type=int,
        nargs="+",
        choices=[0, 4, 6, 8, 12],
        help="the signature for which to compute the upperbound as a space separated list",
    )
    parser.add_argument(
        "--trial_division_bound",
        metavar="bound",
        type=int,
        help="bound up to which to apply trial division for factoring the final result",
        default=10 ** 9,
    )
    parser.add_argument(
        "--aux_bound",
        type=int,
        help="bound on the size of auxiliary primes used",
        default=6,
    )
    parser.add_argument("--verbose", action="store_true", help="get more info printed")
    args = parser.parse_args()
    cli_handler(args)
