# main.py
"""queue_runner.py

Used to run code on all LMFDB number fields.

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


import sys
from sage_code.queue_runner.queue_runner import QueueRunner

if __name__ == "__main__":
    command = sys.argv[1]
    runner = QueueRunner()
    if command == "fill":
        runner.fill_queue(*sys.argv[2:4])
    elif command == "clear":
        runner.clear_queue(sys.argv[2])
    elif command == "len":
        print(runner.len_queue(sys.argv[2]))
    elif command == "run":
        runner.run_function(
            sys.argv[2], "sage_code.lmfdb_isogeny_primes", "isogeny_primes"
        )
    elif command == "dump":
        runner.dump(*sys.argv[2:4])
    elif command == "queue_errors":
        runner.fill_queue_from_errors(sys.argv[2])
