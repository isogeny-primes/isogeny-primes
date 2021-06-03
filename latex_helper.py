"""latex_helper.py

A tool to generate the tables in the introduction of the paper,
minimizing typos.

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait and Maarten Derickx

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

    ====================================================================

"""

import argparse
from sage.all import(Integer, RR, QuadraticField, PolynomialRing, QQ,
                     NumberField, latex)
from isogeny_primes import (get_isogeny_primes, EC_Q_ISOGENY_PRIMES, contains_imaginary_quadratic_field)
from time import perf_counter

import requests  # for accessing LMFDB API
import json

# The URL containing the data we want. The degree will be put into this URL trunk.
URL_TRUNK = ('https://www.lmfdb.org/api/nf_fields/?_format=json&degree={}&'
             'is_galois=true&_fields=label,disc_abs,disc_sign,coeffs')

R = PolynomialRing(QQ, 'x')

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
    the_url = URL_TRUNK.format(str(d))
    data_this_degree = requests.get(url=the_url).json()['data']
    for dat in data_this_degree:
        poly = R(dat['coeffs'])
        K = NumberField(poly, name='a')
        _, contains_hilbert_class_field = contains_imaginary_quadratic_field(K)
        if not contains_hilbert_class_field:
            Delta_K = dat['disc_sign'] * dat['disc_abs']
            return (Delta_K, poly, K)
    raise NotImplementedError("need more hits")

class Latexer(object):

    def __init__(self, degree_range):
        self.range = degree_range

    def lpp_table(self):
        """generate the large putative primes table"""

        output_str = r'${d}$ & ${Delta_K}$ & ${f_K}$ & {possible_isogeny_primes} & ${time_s:.2f}$\\'
        latex_output = []
        for d in range(2, self.range + 1):
            Delta_K, f_K, K = get_smallest_good_number_field(d)
            t_start = perf_counter()
            candidates = get_isogeny_primes(K, norm_bound=50, bound=2000, loop_curves=True)
            t_end = perf_counter()
            time_taken = t_end-t_start
            candidates = [c for c in candidates if c not in EC_Q_ISOGENY_PRIMES]
            candidates.sort()
            possible_isogeny_primes = ", ".join(map(str,candidates))
            f_K_latex = latex(f_K)
            output_here = output_str.format(d=d, Delta_K=Delta_K, f_K=f_K_latex, possible_isogeny_primes=possible_isogeny_primes, time_s=time_taken)
            latex_output.append(output_here)

        for one_line in latex_output:
            print(one_line)


def cli_handler(args):

    DEGREE_RANGE = args.d

    latex_helper = Latexer(DEGREE_RANGE)
    latex_helper.lpp_table()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('d', metavar='d', type=int,
                         help='the d range')
    args = parser.parse_args()
    cli_handler(args)
