/* Genus2Cubic.m

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

function divisors_X_0(N,d)
    assert IsPrime(N);
    tors_order := Numerator((N-1)/12);
    X:=SimplifiedModel(ModularCurveQuotient(N,[]));
    ptsX:=Points(X : Bound:=10000);
    P:=Divisor(ptsX[1]);
    Q:=Divisor(ptsX[2]);
    assert IsPrincipal(tors_order*(P-Q));
    divisors := [];
    n := tors_order div 2;
    for a in [-n + ((tors_order+1) mod 2)..n] do
       D0 := P+Q+(d-2)*Q + a*(P-Q);
       Append(~divisors, D0);
    end for;
    return divisors,[P,Q];
end function;

function HasRealZero(f)
  for fe in Factorization(f) do;
    r,s := Signature(NumberField(fe[1]));
    if r ne 0 then
      return true;
    end if;
  end for;
  return false;
end function;

function IsTotallyNegative(f)
  //don't want to work with corner cases caused by double roots
  assert IsSquarefree(f);
  if HasRealZero(f) then
    return false;
  end if;
  return Evaluate(f,0) lt 0;
end function;

function MySquarefreePart(f);
  //Modifies f by only dividing f by squares
  squarefree_part := SquarefreePart(f);
  square_part := f div squarefree_part;
  f_sq := f div (GCD(square_part,squarefree_part)^2);
  //this assertion fails if a 4th power divides f
  //error out in this edge case
  assert IsSquarefree(f_sq);
  return f_sq;
end function;


function IsELS(C);
  LocalPrimesToCheck := [factorData[1] : factorData in Factorisation(Integers()!Discriminant(C))] cat PrimesUpTo(4*Genus(C)^2);
  for myp in LocalPrimesToCheck do
    if not IsLocallySolvable(C,myp) then
      return false;
    end if;
  end for;

  f,h := HyperellipticPolynomials(SimplifiedModel(C));
  assert h eq 0;
  if IsTotallyNegative(f) then
    return false;
  end if;
  return true;
end function;


S<z> := PolynomialRing(Rationals());

main:=function(PrimesToCheck : twistD := 0);
  badCurves := [];
  for p in PrimesToCheck do
    divisors, cusps := divisors_X_0(p,3);
    R<x,y> := PolynomialRing(Rationals(),2);
    for D in divisors do
      if Degree(Reduction(D,&+cusps)) lt 3 then
        //is of the form hyper elliptic + cusp
        continue;
      end if;
      D1,r := Reduction(D,cusps[1]);
      D1 := D1+r*cusps[1];
      V,f := RiemannRochSpace(D1);
      g := f(V.1);
      assert Degree(g) eq 3;
      eq_g := MinimalPolynomial(g);
      A0 := Parent(Coefficients(g)[1]);
      f := eq_g*(A0 ! LCM([Denominator(c) : c in Coefficients(eq_g)]));
      P := Parent(f);
      P0 := Parent(Coefficient(f,0));
      psi0 := hom< P0 -> R | y>;
      psi := hom< P -> R | psi0,[x]>;
      //psi(f) is now an equation for X_0(p) where g, the map of degree 3 is
      //is just the projection on the x coordinate.
      _,disc := IsUnivariate(Discriminant(psi(f),y));
      f_sq := MySquarefreePart(disc);
      C := HyperellipticCurve(S!f_sq);
      if twistD eq 0 then
        print "genus of C is", Genus(C);
        print p,IsTotallyNegative(f_sq),psi(f);
      else
        twistDSQFree,_ := Squarefree(twistD);
        Ctwist := QuadraticTwist(C, twistDSQFree);
        if IsELS(Ctwist) then
          Append(~badCurves, Ctwist);
          // try to determine the dscriminants. The code in the rest of this block
          // is not currently being used in the main verification for the paper,
          // but could help others doing similar computations.
          rats := RationalPoints(Ctwist : Bound:=50000);
          if #rats gt 0 then
            for my_point in rats do
              theXCoord := my_point[1]/my_point[3];
              psiFX := Evaluate(psi(f),1,theXCoord);
              projZ := hom< R -> S |0,z>;
              defPoly := projZ(psiFX);
              NF<tre>:=NumberField(defPoly);
              theDisc := Discriminant(NF);
              print psiFX, theDisc, Factorisation(theDisc);
            end for;
          end if;
        end if;
      end if;
    end for;
  end for;
  return badCurves;
end function;

//The following computation shows that none of the cubic points on X_0(p) for
//p=23, 29, 31 are totally real, proving Part (1) of Theorem 8.1.

main([23,29,31]);

//The following computation shows that none of the cubic points on X_0(p) for
//p=23, 31 are defined over a cubic field of discriminant equal to -59 up to
// a square. This is required in the proof of Theorem 8.5.

badCurves := main([23,31] : twistD:=-59); // empty, so done.
