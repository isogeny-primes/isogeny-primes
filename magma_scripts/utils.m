/* utils.m
    Some possibly helpful magma functions.

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

function HasCyclicMIsogeny(E,m)

  K := BaseRing(E);
  j := jInvariant(E);
  P<x,y> := PolynomialRing(K,2);
  S<y> := PolynomialRing(K);
  PhiM := P!ClassicalModularPolynomial(m);
  EvalAtj := hom<P -> S | j, y>;

  if HasRoot(EvalAtj(PhiM)) then
    return true;
  else
    return false;
  end if;
end function;

function IsZeroRankofTwist(p,d)

  M := ModularSymbols(p);
  S := CuspidalSubspace(M);
  S_min := AtkinLehnerDecomposition(S)[2];
  chi := BaseExtend(KroneckerCharacter(d),RationalField());
  tw := TwistedWindingElement(S_min,1, chi);
  print tw;
  twmap := RationalMapping(S_min)(tw);
  return not twmap eq Parent(twmap)!0;
end function;
