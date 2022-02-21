/* HigherGenus.m

  Code which does the computations in the final section of the paper.

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

function riemann_roch_spaces_X_0(N,d)
    assert IsPrime(N);
    tors_order := Numerator((N-1)/12);
    X:=SimplifiedModel(ModularCurveQuotient(N,[]));
    ptsX:=Points(X : Bound:=10000);
    P:=Divisor(ptsX[1]);
    Q:=Divisor(ptsX[2]);
    assert IsPrincipal(tors_order*(P-Q));
    dimensions := [];
    unique_divisors := [];
    n := tors_order div 2;
    for a in [-n + ((tors_order+1) mod 2)..n] do
       dimension := [a,Dimension(P+Q+(d-2)*Q + a*(P-Q))];
       Append(~dimensions,dimension);
       if dimension[2] eq 1 then;
           V,phi := RiemannRochSpace(P+Q+(d-2)*Q + a*(P-Q));
           Ds := [D : D in Support(Divisor(phi(V.1))) | Degree(D) eq d];
           if #Ds ne 0 then
               Append(~unique_divisors,Ds[1]);
           end if;
       end if;
    end for;
    return dimensions,unique_divisors;
end function;

function NiceDefiningPolynomial(K)
    O := MaximalOrder(K);
    O := LLL(O);
    fs := [MinimalPolynomial(O.i) : i in [1..Degree(K)] | Degree(MinimalPolynomial(O.i)) eq Degree(K)];
    return fs[1];
end function;

for p in PrimesInInterval(22,72) do
    if p in [37,43,53,61,67] then
        continue;
    end if;
    print p;
    dims,divs := riemann_roch_spaces_X_0(p,2);
    print [dim : dim in dims | dim[2] ne 0];
    dims,divs := riemann_roch_spaces_X_0(p,3);
    fields := [ResidueClassField(D) : D in divs];
    print [dim : dim in dims | dim[2] ne 0],[<Discriminant(K),NiceDefiningPolynomial(K)> : K in fields];
end for;