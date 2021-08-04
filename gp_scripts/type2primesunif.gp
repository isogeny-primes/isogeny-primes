/* partype2primes.gp

    This is the parallelised version of type2primes.gp

    Best to initialise gp like so:

    gp --primelimit 56546719184 --stacksize 20000000

    We gratefully acknowledge the help of Bill Allombert and
    Vinko Petričević in creating this parallelized code. Merci beaucoup
    Bill, and Puno hvala Vinko!

    ====================================================================

    This file is part of Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait and Maarten Derickx

    Isogeny Primes is free software: you can redistribute it and/or
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

\\check if condition CC is satisfied
satisfiesCCunif(p) =
{
  forprime(q = 2,(p/4)^(1/3),
    if(((q^2 + q + 1) % p != 0)&((q^6 + q^3 + 1) % p != 0),
      if(kronecker(q,p) == 1,
          return(0)
        );
      );
    );
  return(1);
}

my_list = List();
blockSize = 100000;

\\print to stdout if p satisfies condition CC
print_satisfiesCCunif(p) =
{
  if(satisfiesCCunif(p),
    listput(my_list,p);
    print(p," is a possible type 2 prime")
  );
}

checktypetwo(pBeg) =
{
    my(p);
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(3,4),
             print_satisfiesCCunif(p));
}

apply(checktypetwo,[0..10]);
listsort(my_list);
print(my_list);
