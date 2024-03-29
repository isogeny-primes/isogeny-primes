/* partype2primes.gp

    This is the parallelised version of type2primes.gp

    Best to initialise gp like so:

    gp --primelimit 56546719184 --stacksize 20000000

    We gratefully acknowledge the help of Bill Allombert and
    Vinko Petričević in creating this parallelized code. Merci beaucoup
    Bill, and Puno hvala Vinko!

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

*/


f = x^3 - x^2 - 2*x + 1;  \\ change this to desired polynomial
K = nfinit(f);
export(K)

typeTwoBound = 1.1e13;  \\ change this to the correct Type 2 bound

\\check if condition CC is satisfied
satisfiesCC(p,q_start) =
{
  forprime(q = q_start,p/4,
    my(odd_fs);
    odd_fs = Set([pqw.f | pqw <- idealprimedec(K,q), pqw.f % 2 == 1]);
    foreach(odd_fs,f,
      if(q^f < p/4,
        if((q^(2*f) + q^f + 1) % p != 0,
          if(kronecker(q,p) == 1,
            return(0);
          );
        );
      );
    );
  );
  return(1);
}
export(satisfiesCC)

\\print to stdout if p satisfies condition CC
print_satisfiesCC(p, q_init) =
{
  if(satisfiesCC(p,q_init),
    print(p," is a type 2 prime")
  );
}
export(print_satisfiesCC)

blockSize=100000;
export(blockSize)


initialCheck(pBeg) =
{
    my(p);
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1, Mod(3,4),
                print_satisfiesCC(p,2));
}
export(initialCheck)

switch=7^(2*poldegree(f)) + 7^(poldegree(f)) + 1;
howManyInitial = floor(switch/blockSize);


my_res_classes = [667,67,547,163,403,43];
export(my_res_classes)

mainCheck(pBeg) =
{
    my(p);
    foreach(my_res_classes,r,
      forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(r,840),
                  print_satisfiesCC(p,11));
    );
}
export(mainCheck)

howMany=floor(typeTwoBound/blockSize);

\\ To check Type 2 primes, run the following two commands

\\ parapply(initialCheck,[0..howManyInitial+1]);
\\ parapply(mainCheck,[howManyInitial..howMany]);
