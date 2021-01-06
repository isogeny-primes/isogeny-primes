# Quadratic Isogeny Primes

This is the repository for the program **Quadratic Isogeny Primes**, whose algorithms are explained in [this paper](https://arxiv.org/).

### What does it do?

You give it a quadratic field, not imaginary quadratic of class number one.

It will give you a set of primes containing the **isogeny primes** for your chosen quadratic field.

It will then be up to you to determine which of those are actually isogeny primes.

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

 - this;

 ```
sage quadratic_isogeny_primes.py --help
```

 - this;

 ```
sage quadratic_isogeny_primes.py --help
```

 - this;

```
sage quadratic_isogeny_primes.py --help
```

### Wasn't there also some Magma and PARI/GP code?

You'll find these in their respective folders.

The [PARI/GP](https://pari-gp.org) code is used to rule out the "type two primes". The file `type2primes.gp` will explain how to do this.

The [Magma](magma.com) code is there so that others can verify the stuff being done in Section 11 of the paper. Note that you'll need a Magma licence for this.

### I found a bug, what do I do now?

Please report any bugs in the [issues](https://github.com/BarinderBanwait/quadratic_isogeny_primes/issues) section.

Alternatively, feel free to send me an email.

### Why only quadratic fields? I want isogeny primes over number fields of degree 29.

Well I'm afraid you're gonna have to wait dear visitor! Even extending the algorithms to cover degree 3 is beyond current technology. One question that needs to be addressed is extending Momose's P_(n) constants; this is related to when the nth Symmetric power of X_0(N) is a formal immersion at infinity. Work that out, and we might be in business!

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

