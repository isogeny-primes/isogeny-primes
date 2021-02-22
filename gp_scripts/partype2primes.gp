/* partype2primes.gp

    This is the parallelised version of type2primes.gp

    Best to initialise gp like so:

    gp --primelimit 56546719184 --stacksize 20000000

    We gratefully acknowledge the help of Bill Allombert and
    Vinko Petričević in creating this parallelized code. Merci beaucoup
    Bill, and Puno hvala Vinko!

    ====================================================================

    This file is part of Quadratic Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait

    Quadratic Isogeny Primes is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    The author can be reached at: barinder.s.banwait@gmail.com

    ====================================================================

*/

\\ 80 billion is the bound for all |D| <= 10

D=-5;  \\ change this to desired value
export(D)

\\check if condition CC is satisfied
satisfiesCC(p) =
{
  forprime(q = 7,p/4,
    if((q^2 + q + 1) % p != 0,
      if(kronecker(D,q) != -1,
        \\ don't memoize the next call
        if(kronecker(q,p) == 1,
          return(0)
        );
      );
    );
  );
  return(1);
}
export(satisfiesCC)

\\print to stdout if p satisfies condition CC
print_satisfiesCC(p) =
{
  if(satisfiesCC(p),
    print(p," is a type 2 prime")
  );
}
export(print_satisfiesCC)

\\for D=3,7
congruence_condition_37(p) =
{
    if(p%24 == 19,
      return(1);
      return(0);
    );
}
export(congruence_condition_37)

\\ for D=-6,-5,6,10
congruence_condition_main(p) =
{
    if(p%24 == 19,
      my(x='x);
      x=p%5;
      if((x==2)||(x==3),
        return(1);
      return(0);
      );
    return(0);
    );
}
export(congruence_condition_main)

\\ for D=5
congruence_condition_5(p) =
{
    if(p%4 == 3,
      x=p%5;
      if((x==2)||(x==3),
        return(1);
      return(0);
      );
    return(0);
    );
}
export(congruence_condition_5)

\\ for D=-10
congruence_condition_m10(p) =
{
    if(p%8 == 3,
      x=p%5;
      if((x==2)||(x==3),
        return(1);
      return(0);
      );
    return(0);
    );
}
export(congruence_condition_m10)

\\for D=2
congruence_condition_2(p) =
{
    if(p%8 == 3,
      return(1);
      return(0);
    );
}
export(congruence_condition_2)

custom_congruence_condition(p,D) =
{
    if (D == -10, congruence_condition_m10(p),
        D == -6, congruence_condition_main(p),
        D == -5, congruence_condition_main(p),
        D == 2, congruence_condition_2(p),
        D == 3, congruence_condition_37(p),
        D == 5, congruence_condition_5(p),
        D == 6, congruence_condition_main(p),
        D == 7, congruence_condition_37(p),
        D == 10, congruence_condition_main(p),
        congruence_condition_main(p));
}
export(custom_congruence_condition)

blockSize=100000;
export(blockSize);

doWork(pBeg) =
{
    my(p,cond);
    forprime(p = pBeg*blockSize, (pBeg+1)*blockSize-1,
             cond=custom_congruence_condition(p,D);
             if(cond,print_satisfiesCC(p)));
}
export(doWork);

howMany=floor(1000000000/blockSize);

\\ Takes about 27 seconds up to a billion!

\\ parapply(doWork,[0..howMany]);
