"""latex_helper.py

A tool to generate the tables in the introduction of the paper,
minimizing typos.

    ====================================================================

    This file is part of Quadratic Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait

    Quadratic Isogeny Primes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    The author can be reached at: barinder.s.banwait@gmail.com

    ====================================================================

"""

import argparse
from sage.all import(Integer, RR, QuadraticField)
from quadratic_isogeny_primes import (CLASS_NUMBER_ONE_DISCS, DLMV,
        get_isogeny_primes, EC_Q_ISOGENY_PRIMES)

class Latexer(object):

    def __init__(self, disc_range):
        self.range = disc_range

    def dlmv_table(self):
        """generate the dlmv table"""

        output_str = r'${Delta_K}$ & $\Q(\sqrt{{{D}}})$ & ${rem} \times 10^{{{exp_at_10}}}$\\'
        for D in range(-self.range,self.range+1):
            if Integer(D).is_squarefree():
                if not D in CLASS_NUMBER_ONE_DISCS:
                    if D != 1:
                        K = QuadraticField(D)
                        Delta_K = K.discriminant()
                        dlmv_bound = RR(DLMV(K))
                        log_dlmv_bound = dlmv_bound.log10()
                        exp_at_10 = int(log_dlmv_bound)
                        rem = log_dlmv_bound - exp_at_10
                        rem = 10**rem
                        rem = rem.numerical_approx(digits=3)
                        output_here = output_str.format(Delta_K=Delta_K, D=D, rem=rem, exp_at_10=exp_at_10)
                        print(output_here)

    def lpp_table(self):
        """generate the large putative primes table"""

        output_str = r'${Delta_K}$ & $\Q(\sqrt{{{D}}})$ & {large_putative_primes}\\'
        latex_output = []
        for D in range(-self.range,self.range+1):
            if Integer(D).is_squarefree():
                if not D in CLASS_NUMBER_ONE_DISCS:
                    if D != 1:
                        K = QuadraticField(D)
                        Delta_K = K.discriminant()
                        candidates = get_isogeny_primes(K, aux_prime_count=4, bound=2000, loop_only_j=True)
                        candidates = [c for c in candidates if c not in EC_Q_ISOGENY_PRIMES]
                        candidates = [c for c in candidates if c > 71]
                        candidates.sort()
                        large_putative_primes = ", ".join(map(str,candidates))
                        output_here = output_str.format(Delta_K=Delta_K, D=D, large_putative_primes=large_putative_primes)
                        latex_output.append(output_here)

        for one_line in latex_output:
            print(one_line)


def cli_handler(args):

    DISC_RANGE = args.R

    latex_helper = Latexer(DISC_RANGE)

    if args.table == 'dlmv':
        latex_helper.dlmv_table()
    elif args.table == 'lpp':
        latex_helper.lpp_table()
    else:
        print("Please specify which table you want to LaTeX, either dlmv or lpp")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('R', metavar='R', type=int,
                         help='the D range')
    parser.add_argument('--table',
                    required=True,
                    nargs='?',
                    choices=['dlmv', 'lpp'],
                    help='which table to generate')

    args = parser.parse_args()
    cli_handler(args)
