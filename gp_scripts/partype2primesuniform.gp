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
  forprime(q = 11,(p/4)^(1/3)+1,
    if(((q^2 + q + 1) % p != 0)&((q^6 + q^3 + 1) % p != 0),
      if(kronecker(q,p) == 1,
          return(0)
        );
      );
    );
  return(1);
}
export(satisfiesCCunif)

\\print to stdout if p satisfies condition CC
print_satisfiesCCunif(p) =
{
  if(satisfiesCCunif(p),
    print(p," is a type 2 prime")
  );
}
export(print_satisfiesCCunif)

blockSize=100000;
export(blockSize);

/*
checktypetwo(pBeg) =
{
    my(p,cond,pm840);
    forprime(p = pBeg*blockSize, (pBeg+1)*blockSize-1,
                 pm840=p%840;
                 cond=(pm840==667)||(pm840==67)||(pm840==547)||(pm840==163)||(pm840==403)||(pm840==43);
                 if(cond,print_satisfiesCCunif(p)));
}
*/

checktypetwo(pBeg) =
{
    my(p);
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(667,840),
                 print_satisfiesCCunif(p));
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(67,840),
                 print_satisfiesCCunif(p));
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(547,840),
                 print_satisfiesCCunif(p));
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(163,840),
                 print_satisfiesCCunif(p));
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(403,840),
                 print_satisfiesCCunif(p));
    forprimestep(p = pBeg*blockSize, (pBeg+1)*blockSize-1,Mod(43,840),
                 print_satisfiesCCunif(p));
}
export(checktypetwo)

\\ howMany=floor(typetwobound/blockSize);

\\ Takes about 5 minutes up to 50 billion! To do that, run the following code:

\\ parapply(checktypetwo,[0/blockSize..50000000000/blockSize]);
