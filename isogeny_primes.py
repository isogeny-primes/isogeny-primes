"""isogeny_primes.py

    Return finite list of isogeny primes attached to a number field.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait and Maarten Derickx

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

    ====================================================================

"""

# Imports

import argparse
import logging

from sage.all import RR, NumberField, exp, log

from sage_code.common_utils import EC_Q_ISOGENY_PRIMES, R
from sage_code.pre_type_one_two import get_pre_type_one_two_primes
from sage_code.type_one_primes import get_type_1_primes
from sage_code.type_two_primes import get_type_2_bound, get_type_2_primes
from sage_code.type_three_not_momose import type_three_not_momose
from sage_code.weeding import apply_weeding

########################################################################
#                                                                      #
#                            DLMV BOUND                                #
#                                                                      #
########################################################################


def DLMV(K):
    """Compute the DLMV bound"""

    # First compute David's C_0

    Delta_K = K.discriminant().abs()
    h_K = K.class_number()
    R_K = K.regulator()
    r_K = K.unit_group().rank()
    delta_K = log(2) / (r_K + 1)
    C_1_K = r_K ** (r_K + 1) * delta_K ** (-(r_K - 1)) / 2
    C_2_K = exp(24 * C_1_K * R_K)
    CHEB_DEN_BOUND = (4 * log(Delta_K ** h_K) + 5 * h_K + 5) ** 2
    C_0 = ((CHEB_DEN_BOUND ** (12 * h_K)) * C_2_K + CHEB_DEN_BOUND ** (6 * h_K)) ** 4

    # Now the Type 1 and 2 bounds

    type_1_bound = (1 + 3 ** (12 * h_K)) ** 2
    type_2_bound = get_type_2_bound(K)

    return max(C_0, type_1_bound, type_2_bound)


########################################################################
#                                                                      #
#                      MAIN CALLING FUNCTION                           #
#                                                                      #
########################################################################


def get_isogeny_primes(
    K,
    norm_bound,
    bound=1000,
    ice_filter=False,
    appendix_bound=1000,
    stop_strategy="auto",
):

    # Start with some helpful user info

    logging.info("Finding isogeny primes for {}.".format(K))
    logging.info("Bound on auxiliary primes is {}.".format(norm_bound))

    # Get and show GenericPrimes

    (strong_type_3_epsilons, embeddings, generic_primes) = get_pre_type_one_two_primes(
        K,
        norm_bound=norm_bound,
        ice_filter=ice_filter,
        stop_strategy=stop_strategy,
    )

    logging.info("generic_primes = {}".format(generic_primes))

    # Get and show TypeOnePrimes

    C_K = K.class_group()

    type_1_primes = get_type_1_primes(K, C_K, norm_bound=norm_bound)
    logging.info("type_1_primes = {}".format(type_1_primes))

    # Get and show TypeTwoPrimes

    type_2_primes = get_type_2_primes(K, embeddings, bound=bound)
    logging.info("type_2_primes = {}".format(type_2_primes))

    # Get and show TypeThreeNotMomosePrimes

    type_3_not_momose_primes, list_of_type_3_fields = type_three_not_momose(
        K, embeddings, strong_type_3_epsilons
    )
    logging.info("type_3_not_momose_primes = {}".format(type_3_not_momose_primes))

    # Put them all together

    candidates = set.union(
        set(type_1_primes),
        set(generic_primes),
        set(type_2_primes),
        set(type_3_not_momose_primes),
    )

    # Try to remove some of these primes via Bruin-Najman and Box tables,
    # Özman sieve, and method of Appendix

    removed_primes = apply_weeding(candidates, K, appendix_bound)

    if removed_primes:
        candidates -= removed_primes
        logging.info("Primes removed via weeding = {}".format(removed_primes))
    else:
        logging.debug("No primes removed via weeding")

    return candidates, list_of_type_3_fields


########################################################################
#                                                                      #
#                            CLI HANDLER                               #
#                                                                      #
########################################################################


def cli_handler(args):  # pylint: disable=redefined-outer-name

    f = R(args.f)

    K = NumberField(f, name="a")

    loglevel = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        level=loglevel,
    )
    logging.debug("Debugging level for log messages set.")

    if args.dlmv:
        dlmv_bound = DLMV(K)
        logging.info(
            "DLMV bound for {} is:\n\n{}\n\nwhich is approximately {}".format(
                K, dlmv_bound, RR(dlmv_bound)
            )
        )
    else:
        logging.warning(
            "Only checking Type 2 primes up to {}. "
            "To check all, use the PARI/GP script.".format(args.bound)
        )
        superset, type_3_fields = get_isogeny_primes(
            K,
            args.norm_bound,
            args.bound,
            args.ice_filter,
            args.appendix_bound,
        )

        superset_list = list(superset)
        superset_list.sort()
        logging.info("superset = {}".format(superset_list))

        possible_new_isog_primes = superset - EC_Q_ISOGENY_PRIMES
        possible_new_isog_primes_list = list(possible_new_isog_primes)
        possible_new_isog_primes_list.sort()
        logging.info(
            "Possible new isogeny primes = {}".format(possible_new_isog_primes_list)
        )
        if type_3_fields:
            how_many_fields = len(type_3_fields)
            logging.info(
                "Outside of the above set, any isogeny primes must be "
                "of Type 3 with imaginary quadratic field L, for L one "
                "of the following {} field(s):\n {}".format(
                    how_many_fields, type_3_fields
                )
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "f", metavar="f", type=str, help="defining polynomial for the Number field"
    )
    parser.add_argument(
        "--norm_bound",
        type=int,
        help="bound on norm of aux primes in PreTypeOneTwo case",
        default=50,
    )
    parser.add_argument("--dlmv", action="store_true", help="get only DLMV bound")
    parser.add_argument(
        "--bound", type=int, help="bound on Type 2 prime search", default=1000
    )
    parser.add_argument(
        "--appendix_bound",
        type=int,
        help="bound on the primes to try the method of the appendix",
        default=1000,
    )
    parser.add_argument("--verbose", action="store_true", help="get more info printed")
    parser.add_argument(
        "--ice_filter",
        action="store_true",
        help="Use the isogeny character enumeration filter",
    )
    args = parser.parse_args()
    cli_handler(args)
