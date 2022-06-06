/* QuadraticVerifs.m

  This verifies the claims in Theorem 1.9/Table 3.

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

/* First run

sage helper_scripts/quadratic.py

at the top level of the repository to obtain the following:

D = -47 possible isogenies = [31, 61]
D = -39 possible isogenies = [97]
D = -37 possible isogenies = [59, 131]
D = -35 possible isogenies = [23]
D = -31 possible isogenies = [73]
D = -23 possible isogenies = [23, 29, 31]
D = -22 possible isogenies = [59]
D = -15 possible isogenies = [23]
D = -5 possible isogenies = [23]
D = 5 possible isogenies = [23, 47]
D = 13 possible isogenies = [31]
D = 17 possible isogenies = [23, 59]
D = 29 possible isogenies = [29]
D = 37 possible isogenies = [23]
D = 41 possible isogenies = [41]
D = 47 possible isogenies = [59]

We will go through each value of D and try to determine whether or not
the list of possible isogenies are isogeny primes or not.

But first some helper functions.

*/

R<x> := PolynomialRing(Rationals());

// simple search on whether there are rational points over K
function HasExtraRatPts(p,D)
    K<a> := NumberField(R![-D,0,1]);
    XK := SmallModularCurve(p,K);
    rats := RationalPoints(XK : Bound:=100);
    return #rats gt 2;
end function;


// If there are no exceptional quadratic points on the modular curve
// then one can try the Twist-Chabauty0 method. If this returns true,
// then there are no non-exceptional QsqrtD points on X0(p). If it
// returns False then it doesn't mean anything.

function TwistChabauty0(p,D)
    X := SmallModularCurve(p);
    X_simp, _ := SimplifiedModel(X);
    Xtwist := QuadraticTwist(X_simp, D);
    Jtwist := Jacobian(Xtwist);
    if RankBound(Jtwist) eq 0 then
        rats := Chabauty0(Jtwist);
        return #rats eq 0;
    else
        // can't conclude, so return
        return false;
    end if;
end function;


// Unfortunately TorsionSubgroup is apparently
// only implemented for Jacobians of genus 2 curves over the rationals,
// so TwistChabauty0 sometimes fails. However we use the method due to
// Samir Siksek mentioned in the arXiv version of `Explicit isogenies of
// prime degree over quadratic fields', wherein if you have a hyperelliptic
// curve whose Jacobian has zero rank and trivial torsion, then the curve
// has no rational points.

function HasTrivialHyperellipticJacobian(p,D)
    X := SmallModularCurve(p);
    X_simp, _ := SimplifiedModel(X);
    Xtwist := QuadraticTwist(X_simp, D);
    Jtwist := Jacobian(Xtwist);
    if IsHyperellipticCurve(Xtwist) then
        if RankBound(Jtwist) eq 0 then
            if TorsionBound(Jtwist,100) eq 1 then
                return true;
            end if;
        end if;
    end if;
    return false;
end function;


/*
    And now for the verifications.
*/


//// D = -47

// 31 determined as follows:

HasExtraRatPts(31,-47); // true


//// D = -31

// This follows from Section 4.7 of Box's paper "Quadratic points on modular
// curves with infinite Mordell-Weil group"


//// D = -23

// 29 and 31 determined as follows:

HasExtraRatPts(29,-23); // true
HasExtraRatPts(29,-23); // true

// 23 determined as follows. No exceptional quadratic points here from
// Bruin-Najman, so can apply TwistChabauty0

TwistChabauty0(23,-23); // true


//// D = -22

// 59 determined as follows.

HasTrivialHyperellipticJacobian(59,-22); // true


//// D = -15

// 23 determined as follows:

HasExtraRatPts(23,-15); // true


//// D = -5

// 23 determined as follows:

HasExtraRatPts(23,-5); // true


//// D = 5

// See the first named author's previous paper.


//// D = 13

// 31 determined as follows:

HasExtraRatPts(31,13); // true


//// D = 17

// 59 determined as follows:

HasTrivialHyperellipticJacobian(59,17); // true


//// D = 29

// 29 determined as follows:

HasExtraRatPts(29,29); // true


//// D = 41

// 41 determined as follows:

HasExtraRatPts(41,41); // true
