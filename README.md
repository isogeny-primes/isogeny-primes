# Quadratic Isogeny Primes

This is the repository for the program **Quadratic Isogeny Primes**, whose algorithms are explained in the paper [Explicit Isogenies of prime degree over Quadratic fields](https://arxiv.org/abs/2101.02673).

### What does it do?

You give it a quadratic field, not imaginary quadratic of class number one.

It will give you a set of primes containing the **isogeny primes** for your chosen quadratic field.

It will then be up to you to determine which of those are actually isogeny primes. Some techniques for how you might go about doing this are explained in Section 11 of the paper.

### How do I use it?

Clone this repo to your computer. It is assumed you have [sage](https://sagemath.org/) installed.

#### Typical use

The main file is `quadratic_isogeny_primes.py`. It takes one positional argument - D - which corresponds to your quadratic field. So if you're interested to see the isogeny primes over Q(root-5), you'd enter the following at the command line:

```
sage quadratic_isogeny_primes.py 5
```

#### Optional arguments

To see the various options, run

```
sage quadratic_isogeny_primes.py --help
```

You'll see that you have the following optional arguments:

 - `--aux_prime_count`; this tells the program how many auxiliary primes to take. So you could for example do the following to take 6 auxiliary primes:

 ```
sage quadratic_isogeny_primes.py 5 --aux_prime_count 6
```

 - `--loop_all`; this will loop over all elliptic curves, not just one for each j-invariant. This is more for sanity checking, since the numbers being computed are twist invariant.

 ```
sage quadratic_isogeny_primes.py 5 --aux_prime_count 6 --loop_all
```

 - `--dlmv`; this will return the dlmv bound:

```
sage quadratic_isogeny_primes.py 17 --dlmv
```

 - `--bound`; this specifies the bound on Type 2 primes that the sage code will check up to. Sage can go up to about 10 million, but beyond that you'll start seeing pari memory overflows.

```
sage quadratic_isogeny_primes.py 17 --bound 10000000
```

 - `--rigorous`; this will attempt to check all primes up to the Type 2 bound. This is a legacy feature which will be phased out in the next release; for checking up to the tens of billions, use the PARI/GP script.

```
sage quadratic_isogeny_primes.py 17 --rigorous
```

#### Running the tests

Some limited unit testing has been implemented. They can be run with the following:

```
sage test_quadratic_isogeny_primes.py
```

Further suggestions on improving the testing, and more test cases, would be very welcome.

#### Generating the tables in the introduction

Typos away!
Your fears we shall allay!
We have code that automatically generates the latex for the tables any day!

The following generates the LaTeX for the DLMV table

```
sage latex_helper.py 10 --table dlmv
```

The `10` refers to how many quadratic fields you want in the table. Change that according to how many you want to see.

Change `dlmv` to `lpp` to get the **Large Putative isogeny Primes** table.

### Wasn't there also some Magma and PARI/GP code?

You'll find these in their respective folders.

The [PARI/GP](https://pari.math.u-bordeaux.fr/) code is used to rule out the "type two primes". The file `type2primes.gp` will explain how to do this.

The [Magma](http://magma.maths.usyd.edu.au/magma/) code is there so that others can verify the stuff being done in Section 11 of the paper. Note that you'll need a Magma licence for this.

### I tried running it for D = 895643215786, but several days later it's still running!

The current release has been tested for |D| <= 100, and even there you'll find cases (D = -86) where the program does not finish in reasonable time. The issue comes, as ever, down to factoring very large numbers.

This is not at all to say that it won't run for |D| > 100; indeed, one of the unit tests checks for D = 2885, which finishes in less than 2 seconds.

Speeding this up to run on massive Ds has been declared a stretch goal for the Birch release.

### I found a bug, what do I do now?

Please report any bugs in the [issues](https://github.com/BarinderBanwait/quadratic_isogeny_primes/issues) section.

If you see anything wrong with specific lines of code, please go to that line, click the three dots that appear, and click "Reference in new issue".

Alternatively, feel free to send me an email.

### Why only quadratic fields? I want isogeny primes over number fields of degree a gazillion.

Well I'm afraid you're gonna have to wait dear reader! Even extending the algorithms to cover degree 3 is beyond current technology. One question that needs to be addressed is extending Momose's P_(n) constants; this is related to when the nth Symmetric power of X_0(N) is a formal immersion at infinity. Work that out, and we might be in business!

### Where can I learn more about elliptic curves and other cool stuff?

You can head over to the [L-functions and Modular Forms Database](https://lmfdb.org/), there you'll find loads of resources, data, and proper cool images.

    ####  Copyright (C) 2021 Barinder Singh Banwait

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
