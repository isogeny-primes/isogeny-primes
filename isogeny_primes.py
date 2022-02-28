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

from sage.all import RR, NumberField, exp, log

from sage_code.common_utils import EC_Q_ISOGENY_PRIMES, R
from sage_code.generic import generic_primes
from sage_code.type_one_primes import type_1_primes
from sage_code.type_two_primes import get_type_2_bound, type_2_primes
from sage_code.type_three_not_momose import type_three_not_momose, type_three_ramified
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
    bound=1000,
    ice_filter=False,
    appendix_bound=1000,
    norm_bound=50,
    auto_stop_strategy=True,
    repeat_bound=4,
):

    # Start with some helpful user info

    logging.info(f"Finding isogeny primes for {K}.")
    logging.info(f"Bound on auxiliary primes is {norm_bound}.")

    # Get and show GenericPrimes

    (strong_type_3_epsilons, embeddings, generic_primes_list) = generic_primes(
        K,
        norm_bound=norm_bound,
        ice_filter=ice_filter,
        auto_stop_strategy=auto_stop_strategy,
        repeat_bound=repeat_bound,
    )

    logging.info(f"generic_primes = {generic_primes_list}")

    # Get and show TypeOnePrimes

    C_K = K.class_group()

    type_1_primes_list = type_1_primes(K, C_K, norm_bound=norm_bound)
    logging.info(f"type_1_primes = {type_1_primes_list}")

    # Get and show TypeTwoPrimes

    type_2_primes_list = type_2_primes(K, embeddings, bound=bound)
    logging.info(f"type_2_primes = {type_2_primes_list}")

    # Get and show TypeThreeNotMomosePrimes

    type_3_not_momose_primes_list, list_of_type_3_fields = type_three_not_momose(
        K, embeddings, strong_type_3_epsilons
    )
    logging.info(f"type_3_not_momose_primes = {type_3_not_momose_primes_list}")

    # Put them all together

    type_three_ramified_list = type_three_ramified(list_of_type_3_fields)

    candidates = set.union(
        set(type_1_primes_list),
        set(generic_primes_list),
        set(type_2_primes_list),
        set(type_3_not_momose_primes_list),
        set(type_three_ramified_list),
    )

    # Try to remove some of these primes via Bruin-Najman and Box tables,
    # Özman sieve, and method of Appendix

    removed_primes = apply_weeding(candidates, K, appendix_bound)

    if removed_primes:
        candidates -= removed_primes
        logging.info(f"Primes removed via weeding = {removed_primes}")
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
            "Only checking Type 2 primes up to %s. "
            "To check all, use the PARI/GP script.",
            args.bound,
        )

        if args.norm_bound:
            norm_bound = args.norm_bound
            auto_stop_strategy = False
        else:
            norm_bound = 50  # still needed for Type 1
            auto_stop_strategy = True

        superset, type_3_fields = get_isogeny_primes(
            K,
            args.bound,
            args.ice,
            args.appendix_bound,
            norm_bound=norm_bound,
            auto_stop_strategy=auto_stop_strategy,
        )

        superset_list = list(superset)
        superset_list.sort()
        logging.info(f"superset = {superset_list}")

        possible_new_isog_primes = superset - EC_Q_ISOGENY_PRIMES
        possible_new_isog_primes_list = list(possible_new_isog_primes)
        possible_new_isog_primes_list.sort()
        logging.info(f"Possible new isogeny primes = {possible_new_isog_primes_list}")
        if type_3_fields:
            how_many_fields = len(type_3_fields)
            logging.info(
                f"Outside of the above set, any isogeny primes must be "
                f"of Momose Type 3 with imaginary quadratic field L, for L one "
                f"of the following {how_many_fields} field(s):\n {type_3_fields}"
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
        "--no_ice",
        dest="ice",
        action="store_false",
        help="Turn off the isogeny character enumeration filter",
    )
    parser.set_defaults(ice=True)
    args = parser.parse_args()
    cli_handler(args)
