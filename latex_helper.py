"""latex_helper.py

A tool to generate the tables in the paper, minimizing typos.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2022 Barinder S. Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    The authors may be reached at: barinder.s.banwait@gmail.com and
    maarten@mderickx.nl.

    ====================================================================

"""

import argparse
import json
from time import perf_counter

import requests  # for accessing LMFDB API
from sage.all import (
    RR,
    NumberField,
    ZZ,
)  # pylint: disable=no-name-in-module
from sage.arith.misc import nth_prime

from isogeny_primes import EC_Q_ISOGENY_PRIMES, get_isogeny_primes, get_type_2_bound

from sage_code.common_utils import CLASS_NUMBER_ONE_DISCS, R

from sage_code.type_one_primes import BAD_FORMAL_IMMERSION_DATA_PATH

# URLs to access LMFDB data
SGNF_URL_TRUNK = (
    "https://www.lmfdb.org/api/nf_fields/?_format=json&degree={}&"
    "is_galois=true&class_number=1&_fields=label,disc_abs,disc_sign,coeffs"
)

LMFDB_NF_URL_TRUNK = "https://www.lmfdb.org/NumberField/{}"


# The following two functions are slightly different to the one in the main code


def _contains_imaginary_quadratic_field_deg_2(K):
    imag_quad = K.is_totally_imaginary()
    hilbert = K.discriminant() in CLASS_NUMBER_ONE_DISCS
    return imag_quad, hilbert


def contains_imaginary_quadratic_field(K):
    """Choosing auxiliary primes in the PreTypeOneTwoCase requires us to
    choose non-principal primes if K contains an imaginary quadratic field."""

    K_deg_abs = K.absolute_degree()

    if K_deg_abs % 2 == 1:
        return False, False

    if K_deg_abs == 2:
        return _contains_imaginary_quadratic_field_deg_2(K)

    quadratic_subfields = K.subfields(2)

    imag_quad_subfields = [
        L for L, _, _ in quadratic_subfields if L.is_totally_imaginary()
    ]

    contains_hilbert_class_field_of_imag_quad = False

    for L in imag_quad_subfields:
        HL = L.hilbert_class_field("c")
        if HL.absolute_degree().divides(K.absolute_degree()):
            K_HL_composite = K.composite_fields(HL)[0]
            if K_HL_composite.absolute_degree() == K_deg_abs:
                contains_hilbert_class_field_of_imag_quad = True
                break

    return (bool(imag_quad_subfields), contains_hilbert_class_field_of_imag_quad)


def get_smallest_good_number_field(d):
    """
    Gets the smallest "good" number field of degree d, where "good" means

    - Galois over Q
    - Does not contain the Hilbert class field of an imaginary quadratic field

    and smallest is in terms of absolute value of discriminant. This connects to
    the LMFDB API for searching, and returns the correct hit.

    Output:

    - Delta_K - discriminant of number field
    - f_K - Definining polynomial of number field

    """
    the_url = SGNF_URL_TRUNK.format(str(d))
    data_this_degree = requests.get(url=the_url).json()["data"]
    for dat in data_this_degree:
        poly = R(dat["coeffs"])
        K = NumberField(poly, name="a")
        Delta_K = dat["disc_sign"] * dat["disc_abs"]
        _, contains_hilbert_class_field = contains_imaginary_quadratic_field(K)
        if not contains_hilbert_class_field:
            return (Delta_K, poly, K, dat["label"])
    raise NotImplementedError("need more hits")


def get_smallest_missing_prime(prime_list):
    """Helper function for the bad formal immersion table"""
    i = 0
    while prime_list[i] == nth_prime(i + 1):
        i += 1
    return nth_prime(i + 1), prime_list[i:]


class Latexer:
    def __init__(self, degree_range):
        self.range = degree_range

    def several_d(self):
        """Table 1 from the introduction"""

        output_str = r"${d}$ & ${Delta_K}$ & {lmfdb_link_latex} & {possible_isogeny_primes} & ${time_s:.2f}$\\"
        latex_output = []
        for d in range(2, self.range + 1):
            Delta_K, _, K, Klabel = get_smallest_good_number_field(d)
            t_start = perf_counter()
            candidates, _ = get_isogeny_primes(K, ice_filter=True, appendix_bound=200)
            t_end = perf_counter()
            time_taken = t_end - t_start
            candidates = [c for c in candidates if c not in EC_Q_ISOGENY_PRIMES]
            candidates.sort()
            possible_isogeny_primes = ", ".join(map(str, candidates))
            lmfdb_link = LMFDB_NF_URL_TRUNK.format(Klabel)
            lmfdb_link_latex = r"\href{{{the_link}}}{{{my_text}}}".format(
                the_link=lmfdb_link, my_text=Klabel
            )
            output_here = output_str.format(
                d=d,
                Delta_K=Delta_K,
                lmfdb_link_latex=lmfdb_link_latex,
                possible_isogeny_primes=possible_isogeny_primes,
                time_s=time_taken,
            )
            print(output_here)
            latex_output.append(output_here)

    def type_2_bounds(self):
        """generate the type 2 bounds table"""

        output_str = (
            r"${d}$ & ${Delta_K}$ & ${label}$ & ${rem} \times 10^{{{exp_at_10}}}$\\"
        )
        latex_output = []
        for d in range(2, self.range + 1):
            Delta_K, f_K, K, label = get_smallest_good_number_field(d)
            type_2_bound = RR(get_type_2_bound(K))
            log_type_2_bound = type_2_bound.log10()
            exp_at_10 = int(log_type_2_bound)
            rem = log_type_2_bound - exp_at_10
            rem = 10**rem
            rem = rem.numerical_approx(digits=3)
            lmfdb_link = LMFDB_NF_URL_TRUNK.format(label)
            lmfdb_link_latex = r"\href{{{the_link}}}{{{my_text}}}".format(
                the_link=lmfdb_link, my_text=label
            )
            output_here = output_str.format(
                d=d,
                Delta_K=Delta_K,
                label=lmfdb_link_latex,
                rem=rem,
                exp_at_10=exp_at_10,
            )
            latex_output.append(output_here)

        for one_line in latex_output:
            print(one_line)

    def bad_formal_immersion_table(self):

        with open(BAD_FORMAL_IMMERSION_DATA_PATH, "r") as bfi_dat_file:
            bfi_dat = json.load(bfi_dat_file)

        degrees = list(bfi_dat.keys())
        degrees_int = [int(d) for d in degrees]
        degrees_int.sort()
        degrees_int = [x for x in degrees_int if x < self.range + 1]
        degrees = [str(d) for d in degrees_int]

        output_str = r"${d}$ & ${smallest_good}$ & {sporadic_bad}\\"
        latex_output = []

        for d in degrees:
            bad_formal_immersion_list = bfi_dat[d]["bad_formal_immersion_list"]

            smallest_good, sporadic_bad_list = get_smallest_missing_prime(
                bad_formal_immersion_list
            )
            sporadic_bad = ", ".join(map(str, sporadic_bad_list))

            output_here = output_str.format(
                d=d, smallest_good=smallest_good, sporadic_bad=sporadic_bad
            )
            latex_output.append(output_here)

        for one_line in latex_output:
            print(one_line)


def bad_formal_immersion_table_old():
    """Not used but still want to keep hold of it for now"""
    with open(BAD_FORMAL_IMMERSION_DATA_PATH, "r") as bfi_dat_file:
        bfi_dat = json.load(bfi_dat_file)

    degrees = list(bfi_dat.keys())
    degrees_int = [int(d) for d in degrees]
    degrees_int.sort()
    degrees = [str(d) for d in degrees_int]

    output_str = r"${d}$ & ${smallest_good}$ & {sporadic_bad} & {almost_good}\\"
    latex_output = []

    for d in degrees:
        bad_formal_immersion_list = bfi_dat[d]["bad_formal_immersion_list"]
        bad_aux_prime_dict = bfi_dat[d]["bad_aux_prime_dict"]

        smallest_good, sporadic_bad_list = get_smallest_missing_prime(
            bad_formal_immersion_list
        )
        sporadic_bad = ", ".join(map(str, sporadic_bad_list))
        bad_aux_prime_dict_rev = {bad_aux_prime_dict[k]: k for k in bad_aux_prime_dict}

        almost_good_primes = list(
            {
                p
                for N in list(bad_aux_prime_dict_rev.keys())
                for p in ZZ(N).prime_divisors()
            }
        )
        almost_good_primes.sort()

        almost_good = {}
        for p in almost_good_primes:
            bad_aux_primes = set()
            for k in bad_aux_prime_dict_rev:
                if k % p == 0:
                    the_bad_aux_prime_here = int(bad_aux_prime_dict_rev[k])
                    bad_aux_primes.add(the_bad_aux_prime_here)
            bad_aux_primes = list(bad_aux_primes)
            bad_aux_primes.sort()
            bad_aux_primes_str = ", ".join(map(str, bad_aux_primes))
            almost_good[p] = bad_aux_primes_str

        almost_good = ", ".join(
            ["{} ({})".format(p, q) for p, q in almost_good.items()]
        )

        output_here = output_str.format(
            d=d,
            smallest_good=smallest_good,
            sporadic_bad=sporadic_bad,
            almost_good=almost_good,
        )
        latex_output.append(output_here)

    for one_line in latex_output:
        print(one_line)


def cli_handler(args):

    THE_RANGE = args.d

    latex_helper = Latexer(THE_RANGE)

    if args.table == "several_d":
        latex_helper.several_d()
    elif args.table == "lpp":
        latex_helper.lpp_table()
    elif args.table == "bfi":
        latex_helper.bad_formal_immersion_table()
    elif args.table == "t2b":
        latex_helper.type_2_bounds()
    else:
        print("Please specify which table you want to LaTeX")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("d", metavar="f", type=int, help="the range depending on table")
    parser.add_argument(
        "--table",
        required=True,
        nargs="?",
        choices=["several_d", "lpp", "bfi", "t2b"],
        help="which table to generate",
    )
    args = parser.parse_args()
    cli_handler(args)
