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

from sage.all import Integer, prime_divisors

from sage_code.strong_uniform_bounds import unif_bd, type_one_unif_primes


def do_uniform_quadratic():

    EPSILONS = [(0,4), (0,8), (4,4), (4,6), (4,12)]

    for eps in EPSILONS:
        bd = Integer(unif_bd(2,eps))
        logging.info(f"Isogeny primes for {eps} are {bd.prime_divisors()}")

    type_one_bound = type_one_unif_primes(int(2))

    logging.info(f"Type 1 uniform bound is {type_one_bound}")


def do_uniform(d, eps):
    tr_not_6_bd = Integer(unif_bd(d, eps))
    logging.info(f"Trace not 6 mod 12 bound is {tr_not_6_bd}")
    # logging.info(f"Trace not 6 mod 12 primes are {tr_not_6_bd.prime_divisors()}")

    type_one_bound = type_one_unif_primes(d)

    logging.info(f"Type 1 uniform bound is {type_one_bound}")


def get_eps_from_input(tuple_as_str):
    return tuple([Integer(t) for t in tuple_as_str.split(",")])


def cli_handler(args):  # pylint: disable=redefined-outer-name

    loglevel = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        level=loglevel,
    )
    logging.debug("Debugging level for log messages set.")

    eps = get_eps_from_input(args.eps)

    if args.d == int(2):
        do_uniform_quadratic()
    else:
        do_uniform(args.d, eps)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "d",
        metavar="d",
        type=int,
        help="bound on the primes to try the method of the appendix",
    )
    parser.add_argument(
        "--eps",
        required=True,
        type=str,
        help="bound on the primes to try the method of the appendix",
    )
    parser.add_argument("--verbose", action="store_true", help="get more info printed")
    args = parser.parse_args()
    cli_handler(args)