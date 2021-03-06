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
from isogeny_primes import (get_isogeny_primes, EC_Q_ISOGENY_PRIMES,
    contains_imaginary_quadratic_field, CLASS_NUMBER_ONE_DISCS)
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


# from isogeny_primes import (get_isogeny_primes, EC_Q_ISOGENY_PRIMES,
#     contains_imaginary_quadratic_field, CLASS_NUMBER_ONE_DISCS)
# URL_TRUNK = ('https://www.lmfdb.org/api/nf_fields/?_format=json&degree={}&'
#              'is_galois=false&_fields=label,disc_abs,disc_sign,coeffs')
# norm_bound = 50
# R = PolynomialRing(QQ, 'x')
# def find_stop_iteration(d):
#     the_url = URL_TRUNK.format(str(d))
#     data_this_degree = requests.get(url=the_url).json()['data']
#     for dat in data_this_degree:
#         poly = R(dat['coeffs'])
#         K = NumberField(poly, name='a')

#         Kgal = K.galois_closure('b')
#         aux_primes =[]
#         contains_imaginary_quadratic, contains_hilbert_class_field = contains_imaginary_quadratic_field(Kgal)

#         try:
#             it = Kgal.primes_of_degree_one_iter(max_iterations=1000)
#             aux_prime_count = 2
#             while aux_prime_count > 0:
#                 aux_prime_candidate = next(it)
#                 if (not contains_imaginary_quadratic) or (not aux_prime_candidate.is_principal()):
#                     if aux_prime_candidate.norm() > norm_bound:
#                         aux_primes.append(aux_prime_candidate)
#                     aux_prime_count -= 1
#         except StopIteration:
#             print(K, dat['label'])


class Latexer(object):

    def __init__(self, degree_range):
        self.range = degree_range

    def several_d(self):
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

    def lpp_table(self):
        """generate the large putative primes table"""

        output_str = r'${Delta_K}$ & $\Q(\sqrt{{{D}}})$ & {new_isogeny_primes}\\'
        latex_output = []
        for D in range(-self.range,self.range+1):
            if Integer(D).is_squarefree():
                if not D in CLASS_NUMBER_ONE_DISCS:
                    if D != 1:
                        K = QuadraticField(D)
                        Delta_K = K.discriminant()
                        candidates = get_isogeny_primes(K, norm_bound=150, bound=2000, loop_curves=True)
                        candidates = [c for c in candidates if c not in EC_Q_ISOGENY_PRIMES]
                        # candidates = [c for c in candidates if c > 71]
                        candidates.sort()
                        new_isogeny_primes = ", ".join(map(str,candidates))
                        output_here = output_str.format(Delta_K=Delta_K, D=D, new_isogeny_primes=new_isogeny_primes)
                        latex_output.append(output_here)
                        print(output_here)
        print("Total number of hits: {}".format(len(latex_output)))
        # for one_line in latex_output:
        #     print(one_line)

def cli_handler(args):

    THE_RANGE = args.d

    latex_helper = Latexer(THE_RANGE)

    if args.table == 'several_d':
        latex_helper.several_d()
    elif args.table == 'lpp':
        latex_helper.lpp_table()
    else:
        print("Please specify which table you want to LaTeX, either dlmv or lpp")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('d', metavar='d', type=int,
                         help='the range depending on table')
    parser.add_argument('--table',
                required=True,
                nargs='?',
                choices=['several_d', 'lpp'],
                help='which table to generate')
    args = parser.parse_args()
    cli_handler(args)
