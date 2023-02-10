# Isogeny Primes

This is the repository for the program **Isogeny Primes** explained in the papers [Explicit isogenies of prime degree over number fields](https://arxiv.org/abs/2203.06009) and `Towards strong uniformity for isogenies of prime degree'.

### What does it do?

There are two ways of using the program. 

1. In **Uniform** mode. This takes an input polynomial, and will give you a set of primes containing the **isogeny primes** for the number field defined by your input polynomial.

2. In **Strong Uniform** mode. This takes an input integer, and will give you a set of primes containing the **isogeny primes** for *any number field of degree d* broken down by the trace of the signature of the isogeny. Not all trace values can be dealt with currently; see the `Towards strong uniformity' paper for more on what can currently be handled.


It will then be up to you to determine which of those output lists are actually isogeny primes.

### How do I use it?

Clone this repo to your computer. It is assumed you have [sage](https://sagemath.org/) installed.

#### Typical use - Uniformity

The main file is `isogeny_primes.py`. It takes one positional argument - f - which is the polynomial. So if you're interested to see the isogeny primes over Q(zeta7), you'd enter the following at the command line:

```
sage isogeny_primes.py 'x^6 - x^5 + x^4 - x^3 + x^2 - x + 1'
```

#### Typical use - Strong Uniformity

The main file is `uniform_isogeny_primes.py`. It takes one positional argument - d - which is the degree. So if you're interested to see the isogeny primes over cubic fields, you'd enter the following at the command line:

```
sage uniform_isogeny_primes.py 3
```

#### That didn't work for me ...

**Isogeny Primes** has a number of requirements which your system may not have. In this case, you should run `make` in the top level of the directory. This will install the necessary requirements. You could also do `make valid_fast`; this will install everything _and_ run some light testing to make sure everything works fine. If anything goes wrong and you want to start again, run `make clean` before trying again.

#### Optional arguments

To see the various options, run

```
sage isogeny_primes.py --help
```

or ```
sage uniform_isogeny_primes.py --help
```

For the first command you'll see that you have the following optional arguments:

 - `--norm_bound`; if specified this will take auxiliary primes up to the specified bound. However, doing this also turns off the auto stop strategy.

 ```
sage isogeny_primes.py 'x^3 - 5' --norm_bound 50
```

 - `--dlmv`; returns the DLMV bound for the number field.

```
sage isogeny_primes.py 'x^3 - x^2 + 5*x + 14' --dlmv
```

 - `--bound`; this specifies the bound on Type 2 primes that the sage code will check up to. Sage can go up to about 10 million, but beyond that you'll start seeing pari memory overflows.

```
sage isogeny_primes.py 'x^5 + 19' --bound 10000000
```

 - `--appendix_bound`; the bound on possible isogeny primes beyond which the method of the appendix will not be performed. Increasing this can massively slow down the code.

```
sage isogeny_primes.py 'x^6 + 6*x + 14*x + 89' --appendix_bound 1000
```

 - `--verbose`; for power users who want to know what the program is doing.

```
sage isogeny_primes.py 'x^5 + 19' --verbose
```

 - `--no_ice`; this turns off the isogeny character enumeration (ICE) filter. This will speed up the code but likely give a larger superset. Turning ICE off for large degree number fields might be necessary.

```
sage isogeny_primes.py 'x^17 + 19' --no_ice
```

#### Generating the tables in the paper

Code for generating the LaTeX for the tables in the paper is done with `latex_helper.py`. For example,

```
sage latex_helper.py 10 --table bfi
```

will generate the Bad Formal Immersion Data table. The `10` refers to how many fields you want in the table. Change that according to how many you want to see.

See the source code of `latex_helper.py` to see the other options.


### I found a bug, what do I do now?

Please report any bugs in the [issues](https://github.com/isogeny_primes/isogeny_primes/issues) section.

If you see anything wrong with specific lines of code, please go to that line, click the three dots that appear, and click "Reference in new issue".

Alternatively, feel free to send us an email.

# Developer's Guide

## Project layout
The directory layout is as follows

    .
    ├── gp_scripts          # scripts written in pari/gp
    ├── helper_scripts      # scripts helping with starting sage
    ├── magma_scripts       # magma scripts
    ├── sage_code           # main python/sage source directory
    │   └── data_files      # json files containing static data and caches
    └── tests
        ├── fast_tests      # unittests - these test should be fast and match
        │                   #             the sage_code file structure
        └── slow_tests      # integrationtests - these tests take longer


## Developer tools

The most important developer tool is the Makefile. It contains several commands that help with development.

The most important command are:

    make valid
    make valid_fast

At least one of these should be run before pushing your code to github and
creating a pull request. Both of these commands do the following things

- Make sure all runtime and development dependencies are installed
- Fix all code style issues using black
- Run several code quality inspection tools
- Run (a subset of) tests in the tests folder
- Print a coverage report of tests that were run

The main difference is that `make valid` runs all tests, while `make valid_fast` runs only unittests.

Here are some other useful commands:

Create a virtual environment that uses the python in sage as base intepreter:

    make venv

Install all development requirements as specified in the requirements.in
and the requirements-dev.in files:

    make pip-install-dev

Run all tests:

    make tests

For the rest see the source code of the makefile itself.

# Copyright

    ####  Copyright (C) 2022 Barinder S. Banwait and Maarten Derickx

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