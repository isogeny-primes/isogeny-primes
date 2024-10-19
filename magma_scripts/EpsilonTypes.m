/*
The main function in this file is eps_from_isogeny. Which computes epsilon types of isogenies.
It is inspired by section 1 of Serre's `Proprietes Galoisienne' paper.

To generate the curves in Proposition 7.3 of the paper, run the following:

load "EpsilonTypes.m";
print_eps_type_info(0,17);
print_eps_type_info(0,19);
print_eps_type_info(2,1: special_mode:="37");

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



function eps_from_isogeny_raw(phi, p);
  phi_formal := FormalGroupHomomorphism(phi,3);
  e := Valuation(Coefficient(phi_formal,1),p);
  return e;
end function;

function eps_from_isogeny_local(phi, pp);
  p := Degree(phi);
  assert Norm(pp) mod p eq 0;

  E := Domain(phi);
  E2 := Codomain(phi);
  j := jInvariant(E);
  vp_j := Valuation(j,pp);
  OK := Order(pp);
  K := NumberField(OK);
  R<U> := PolynomialRing(K);

  //i haven't though about the case where this doesn't hold
  assert RamificationIndex(pp) eq 1;


  exp := 2;
  if vp_j gt 0 then
    exp := 6;
  end if;
  if Valuation(j-1728,pp) gt 0 then
    exp := 4;
  end if;
  if vp_j lt 0 then
    j2 := jInvariant(E2);
    vp_j2 := Valuation(j2,pp);
    // this case corresponds to K^*/q^pZ -> K^*/(zeta_p^i*q)^Z
    if  vp_j/vp_j2 eq p then return 0; end if;
    // this case corresponds to K^*/q^Z -> K^*/<zeta_p,q^Z> ~= K^*/q^(pZ)
    if  vp_j/vp_j2 eq 1/p then return 12; end if;
    print vp_j, vp_j2, j, j2;
    error Error("Potential multiplicative reduction went wrong valuation ratio should be p or 1/p!");
  end if;



  L<u> := NumberField(U^exp-p);
  OL := MaximalOrder(L);
  fac := Factorization(ideal< OL | [OL ! x : x in Generators(pp)] >);
  ppL := fac[1][1];
  assert #fac eq 1;
  assert fac[1][2] eq exp;
  EL := BaseChange(E,L);
  EL2 := BaseChange(E2,L);
  locL, ELp := LocalInformation(EL,ppL);
  f := IsomorphismToIsogeny(Isomorphism(ELp,EL));
  loc2L, EL2p := LocalInformation(EL2,ppL);
  f2 := IsomorphismToIsogeny(Isomorphism(EL2,EL2p));
  e1 := eps_from_isogeny_raw(f, ppL);
  e2 := exp*eps_from_isogeny_raw(phi, pp);
  e3 := eps_from_isogeny_raw(f2, ppL);
  //print e1,e2,e3;
  return (12 div exp)*(e1 + e2 + e3);

end function;

function eps_from_isogeny(phi);
  /*
  The input should be an isogeny of prime degree p, over a numberfield K.
  Currently the code only works if p is unramified in K.
  It might be possible to extend it to the case where p is larger then
  its ramification index in K.
  */
  p := Degree(phi);
  assert IsPrime(p);
  E := Domain(phi);
  K := BaseField(E);
  OK := MaximalOrder(K);
  eps_type := &cat[
    [eps_from_isogeny_local(phi,p_e[1]) : i in [1..InertiaDegree(p_e[1])]]
    : p_e in Factorization(p*OK)
  ];
  return eps_type;
end function;


/*
from this point on the code is just for generating isogenies and printing info about them
*/

function quadratic_isogeny(t,p)
  /*
  This function generates quadratic points on the modular curve X_0(p), when p is such that this curve is (hyper) elliptic.
  It needs a rational number t as input. The quadratic point will have t as image under the (hyper) elliptic map.

  The output is a tuple consisting of the isogeny and the maximal order in the quadratic field.
  */
  R<Y> := PolynomialRing(RationalField());
  X0 := SmallModularCurve(p);
  assert (Genus(X0) eq 1) or IsHyperelliptic(X0);
  j := jInvariant(X0,p);
  if Genus(X0) eq 1 then
    C,phi := WeierstrassModel(X0);
  else
    C,phi := SimplifiedModel(X0);
  end if;
  f := HyperellipticPolynomials(C);
  y2 := Evaluate(f,t);
  K<y> := NumberField(Y^2-y2);
  O_K := MaximalOrder(K);
  P_C :=  C(K) ! [t,y];
  P := (phi^(-1))(P_C);
  phi := Isogeny(P,p);
  return phi,O_K;
end function;

function quadratic_isogeny_37(a, b)
  /*
  This function generates quadratic points on the modular curve X_0(37). The mordell weil group of J_0(37) is generated
  by two elements x and y. So any point of Pic^{2} X_0(37) is of the form  K + a*x + b*y, where a and b are integers and
  K is a canonical divisor. Outside the hyperelliptic locus the map X_0(37)^{(2)}(Q) -> Pic^{2} X_0(37)(Q) is a
  bijection.  The tuple (a,b) should not be (0,0). The isogeny returned is the one corresponding to unique point of
  X_0(37)^{(2)}(Q) mapping to K + a*x + b*y in Pic^{2} X_0(37)(Q).

  The output is a tuple consisting of the isogeny and the maximal order in the quadratic field.
  */
  assert [a mod 3, b] ne [0 mod 3, 0];
  X0 := SmallModularCurve(37);
  J0 := Jacobian(X0);
  //we should probably cache the output of the following command;
  MW,f,b1,b2 := MordellWeilGroup(J0);
  assert b1;
  assert b2;
  P := Eltseq(f(a*MW.1+b*MW.2));
  assert Degree(P[1]) eq 2;
  assert IsIrreducible(P[1]);
  K := NumberField(P[1]);
  x0 := K.1;
  y0 := Evaluate(P[2], x0);
  P0 := X0(K) ! [x0, y0];
  phi := Isogeny(P0, 37);
  O_K := MaximalOrder(K);
  return phi,O_K;
end function;


function print_eps_type_info(t, p: special_mode:="None")
  if special_mode eq "None" then
    phi,OK := quadratic_isogeny(t,p);
  elif special_mode eq "37" then
    phi,OK := quadratic_isogeny_37(t,p);
  else
    print("special_mode should be 'None' or '37'");
    assert false;
  end if;

  K := NumberField(OK);
  print t, ClassNumber(K), K;
  E := Domain(phi);
  E2 := Codomain(phi);
  j := jInvariant(E);
  j2 := jInvariant(E2);
  D_E := Discriminant(E);
  print "found isogeny between j1 =", j, "E1 =", E;
  print "and j2 =", j2, "E2 =", E2;
  for pp_e in Factorization(p*OK) do;
    pp := pp_e[1];
    loc,Ep := LocalInformation(E,pp);
    Dp := Discriminant(Ep);
    print t,pp_e[2],Valuation(j,pp),Valuation(j-1728,pp),Valuation(Dp,pp),loc[5];


  end for;
  eps := eps_from_isogeny(phi);
  print "epsilon type:", eps;
  return j,j2, eps;
end function;

/*
p := 11;
for t in [3,7,10] do
  tmp := print_eps_type_info(t, p);
end for;

p := 23;
for t in [3,7,10,17,20] do
  tmp := print_eps_type_info(t, p);
end for;
*/
