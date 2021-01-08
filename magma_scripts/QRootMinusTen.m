/* QRootMinusTen.m
    Magma code related to the final section of the paper, for Q(\sqrt{-10})

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
K<a> := NumberField(R![10, 0, 1]);

// The only cases done with Magma are 29 and 31

// 29 - Bruin-Najman-Twist-Chabauty0 Method
N:=29;
X := SmallModularCurve(N,K);
rats := RationalPoints(X : Bound:=300);
rats;  // trivial, so change tack
X := SmallModularCurve(N);
X_simp, f := SimplifiedModel(X);
X7 := QuadraticTwist(X_simp,-10);
RationalPoints(X7 : Bound:=100000);
J7 := Jacobian(X7);
RankBound(J7);
Chabauty0(J7);

// 31 - Bruin-Najman-Twist-Chabauty0 Method
N:=31;
X := SmallModularCurve(N,K);
rats := RationalPoints(X : Bound:=300);
rats;  // trivial, so change tack
X := SmallModularCurve(N);
X_simp, f := SimplifiedModel(X);
X7 := QuadraticTwist(X_simp,-10);
RationalPoints(X7 : Bound:=100000);
J7 := Jacobian(X7);
RankBound(J7);
Chabauty0(J7);
