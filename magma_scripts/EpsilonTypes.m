/*
The main function in this file is eps_from_isogeny. Which computes epsilon types of isogenies.
It is inspired by section 1 of:
    https://www.wstein.org/sage_summer/bsd_comp/Serre-properties_galoisiennes_des_points_dordre_fini_des_courbes_elliptiques.pdf
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
  OK := Order(pp);
  K := NumberField(OK);
  R<U> := PolynomialRing(K);

  exp := 2;
  if Valuation(j,pp) gt 0 then
    exp := 6;
  end if;
  if Valuation(j-1728,pp) gt 0 then
    exp := 4;
  end if;
  if Valuation(j,pp) lt 0 then
    error Error("Potential multiplicative reduction not implemented yet!");
  end if;

  //so we can just take roots out of p
  //in order to get semistable reduction
  assert RamificationIndex(pp) eq 1;

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

  As output it just gives the j-invariant of this quadratic point.
  */
  R<Y> := PolynomialRing(RationalField());
  x0 := SmallModularCurve(p);
  assert (Genus(x0) eq 1) or IsHyperelliptic(x0);
  j := jInvariant(x0,p);
  if Genus(x0) eq 1 then
    C,phi := WeierstrassModel(x0);
  else
    C,phi := SimplifiedModel(x0);
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



function print_eps_type_info(t, p)
  phi,OK := quadratic_isogeny(t,p);
  K := NumberField(OK);
  print t,K;
  E := Domain(phi);
  j := jInvariant(E);
  j2 := jInvariant(Codomain(phi));
  D_E := Discriminant(E);
  print "found isogeny between j1 =", j;
  print "and j2 =", j2;
  for pp_e in Factorization(p*OK) do;
    pp := pp_e[1];
    loc,Ep := LocalInformation(E,pp);
    Dp := Discriminant(Ep);
    print t,t mod p,pp_e[2],Valuation(j,pp),Valuation(j-1728,pp),Valuation(Dp,pp),loc[5];


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
