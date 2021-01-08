/* QRootMinusFive.m
    Magma code related to the final section of the paper, for Q(\sqrt{-5})

    ====================================================================

    This file is part of Quadratic Isogeny Primes.

    Copyright (C) 2021 Barinder Singh Banwait

    Quadratic Isogeny Primes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
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

// Globals
R<x> := PolynomialRing(Rationals());
K<a> := NumberField(R![5, 0, 1]);

// 23 - Magma Point Search Method
N:=23;
X := SmallModularCurve(N,K);
rats := RationalPoints(X : Bound:=10);
P := rats[3];
PConj := rats[4];
jInvariant(P,N);
jInvariant(PConj,23);

// 29 - Özman sieve - done with sage code

// 31 - Bruin-Najman-Twist-Chabauty0 Method
N:=31;
X := SmallModularCurve(N,K);
rats := RationalPoints(X : Bound:=300);
rats;  // trivial, so change tack
X := SmallModularCurve(N);
X_simp, f := SimplifiedModel(X);
X5 := QuadraticTwist(X_simp,-5);
RationalPoints(X5 : Bound:=50000);
J5 := Jacobian(X5);
RankBound(J5);
Chabauty0(J5);

// 41, 47, 53, 59, 61, 71, 73 do not use Magma code