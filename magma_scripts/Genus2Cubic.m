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
        Append(~badCurves, C);
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
              if my_point[3] eq 0 then continue; end if;
              theXCoord := my_point[1]/my_point[3];
              psiFX := Evaluate(psi(f),1,theXCoord);
              projZ := hom< R -> S |0,z>;
              defPoly := projZ(psiFX);
              if Degree(defPoly) ne 3 or not IsIrreducible(defPoly) then continue; end if;
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


/* analyze_squareclass(p, twistD)
 *
 * Per-squareclass analysis for prime p and signed (squarefree-or-not) twistD.
 * Returns < n_els, witnesses > where:
 *   n_els     :: Integer  -- number of ELS-surviving discriminant-curve twists
 *   witnesses :: List of tuples < poly_coeffs, theDisc, fac_str >
 *                  poly_coeffs :: list of integers (defining poly, ascending,
 *                                 denominators cleared)
 *                  theDisc     :: Integer (disc of the recovered cubic field)
 *                  fac_str     :: String (factorisation of theDisc)
 *
 * Empty `witnesses` with `n_els = 0` means provable NO (every relevant twist
 * fails ELS) for any field with disc(K) = twistD up to squares.
 */
analyze_squareclass:=function(p, twistD)
  n_els := 0;
  witnesses := [* *];

  divisors, cusps := divisors_X_0(p, 3);
  R<x,y> := PolynomialRing(Rationals(), 2);

  twistDSQFree, _ := Squarefree(twistD);

  for D in divisors do
    if Degree(Reduction(D, &+cusps)) lt 3 then
      continue;
    end if;
    D1, r := Reduction(D, cusps[1]);
    D1 := D1 + r * cusps[1];
    V, fmap := RiemannRochSpace(D1);
    g := fmap(V.1);
    assert Degree(g) eq 3;
    eq_g := MinimalPolynomial(g);
    A0 := Parent(Coefficients(g)[1]);
    f := eq_g * (A0 ! LCM([Denominator(c) : c in Coefficients(eq_g)]));
    P := Parent(f);
    P0 := Parent(Coefficient(f, 0));
    psi0 := hom< P0 -> R | y >;
    psi := hom< P -> R | psi0, [x] >;
    _, disc := IsUnivariate(Discriminant(psi(f), y));
    f_sq := MySquarefreePart(disc);
    C := HyperellipticCurve(S ! f_sq);

    Ctwist := QuadraticTwist(C, twistDSQFree);
    if not IsELS(Ctwist) then
      continue;
    end if;
    n_els +:= 1;

    rats := RationalPoints(Ctwist : Bound := 50000);
    for my_point in rats do
      if my_point[3] eq 0 then continue; end if;
      theXCoord := my_point[1] / my_point[3];
      psiFX := Evaluate(psi(f), 1, theXCoord);
      projZ := hom< R -> S | 0, z >;
      defPoly := projZ(psiFX);
      if Degree(defPoly) ne 3 or not IsIrreducible(defPoly) then continue; end if;
      denom_lcm := LCM([Denominator(c) : c in Coefficients(defPoly)]);
      defPolyZ := defPoly * denom_lcm;
      poly_coeffs := [Integers() ! c : c in Coefficients(defPolyZ)];
      NF<tre> := NumberField(defPoly);
      theDisc := Discriminant(NF);
      fac_str := Sprintf("%o", Factorisation(theDisc));
      Append(~witnesses, < poly_coeffs, theDisc, fac_str >);
    end for;
  end for;

  return n_els, witnesses;
end function;


/* run_pipeline(input_csv, output_csv, prime, lim)
 *
 * Reads input_csv (header "D"; one signed integer per row), runs
 * analyze_squareclass(prime, D) per row, and writes output_csv with header:
 *   D,n_els_survivors,n_witnesses,witness_disc_list,witness_poly_list,runtime_s
 *
 * Encoding:
 *   witness_disc_list :: ";"-separated signed ints, e.g. "-107;-428"
 *   witness_poly_list :: ";"-separated coeffs (ascending) within a poly,
 *                        "|"-separated across polys, e.g. "32;-48;16;16|4;-12;4;4"
 *
 * lim > 0 truncates to the first lim rows (smoke-test convenience).
 */
run_pipeline:=procedure(input_csv, output_csv, prime, lim)
  SetColumns(0);

  printf "Input  : %o\n", input_csv;
  printf "Output : %o\n", output_csv;
  printf "Prime  : %o\n", prime;
  if lim gt 0 then
    printf "Limit  : %o\n", lim;
  end if;
  printf "\n";

  input_lines := Split(Read(input_csv), "\n");
  input_rows := [l : l in input_lines[2..#input_lines] | l ne ""];
  if lim gt 0 and lim lt #input_rows then
    input_rows := input_rows[1..lim];
  end if;

  printf "Processing %o squareclasses.\n\n", #input_rows;

  PrintFile(output_csv,
    "D,n_els_survivors,n_witnesses,witness_disc_list,witness_poly_list,runtime_s"
    : Overwrite := true);

  t_total := Cputime();

  for i := 1 to #input_rows do
    fields := Split(input_rows[i], ",");
    D := StringToInteger(fields[1]);

    try
      t0 := Cputime();
      n_els, witnesses := analyze_squareclass(prime, D);
      dt := Cputime(t0);

      disc_list_str := Join([Sprintf("%o", w[2]) : w in witnesses], ";");
      poly_list_str := Join([Join([Sprintf("%o", c) : c in w[1]], ";")
                             : w in witnesses], "|");

      row := Sprintf("%o,%o,%o,%o,%o,%.2o",
        D, n_els, #witnesses, disc_list_str, poly_list_str, dt);
      PrintFile(output_csv, row);

      printf "  [%3o/%o] D=%-7o n_els=%o n_wit=%o (%.2os)\n",
        i, #input_rows, D, n_els, #witnesses, dt;
    catch e
      err_msg := Sprintf("%o", e);
      printf "  [%3o/%o] D=%-7o ERROR: %o\n", i, #input_rows, D, err_msg;
      PrintFile(output_csv, Sprintf("%o,,,,\"ERROR: %o\",", D, err_msg));
    end try;
  end for;

  printf "\nTotal time: %.2o s\n", Cputime(t_total);
end procedure;


//The following computation shows that none of the cubic points on X_0(p) for
//p=23, 29, 31 are totally real, proving Part (1) of Theorem 9.1.
//
//The following computation shows that none of the cubic points on X_0(p) for
//p=23, 31 are defined over a cubic field of discriminant equal to -59 up to
// a square. This is required in the proof of Theorem 9.5.
//
//Both paper-verification calls are guarded so this file can be loaded as a
//library without auto-running them. To rerun the paper checks, invoke as:
//  magma RUN_PAPER_CHECKS:=1 magma_scripts/Genus2Cubic.m

if assigned RUN_PAPER_CHECKS then
  main([23,29,31]);
  badCurves := main([23,31] : twistD:=-59); // empty, so done.
end if;
