#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]

"Here, we specialize the script 'TQF ring calcs' to the case where we have imposed various divisibility/coprimality conditions.
We write r=8r', u=2u', b=4b',c=16c', d=64d', where 4det(Ax+By) has coeffs (a,b,c,d)."
"Basis of Z[alpha] given by powers of alpha"

m0 = IdentityMatrix[4];
m1 = {{0,0,0,(q^2-(b^2-3a c)y^2)/12},{1,0,0,-r},{0,1,0,-q},{0,0,1,0}};
m2 = {{0,0,(q^2-(b^2-3a c)y^2)/12,0},{0,0,-r,(q^2-(b^2-3a c)y^2)/12},{1,0,-q,-r},{0,1,0,-q}};
m3 = {{0,(q^2-(b^2-3a c)y^2)/12,0,-q (q^2-(b^2-3a c)y^2)/12},{0,-r,(q^2-(b^2-3a c)y^2)/12,q r},{0,-q,-r,(13 q^2-(b^2-3a c)y^2)/12},{1,0,-q,-r}};

m = m3 + A m2 + B m1 + F m0;
m//MatrixForm;
p = CharacteristicPolynomial[m,w]//Expand;
Coefficient[p,w^4];
Coefficient[p,w^3];
Coefficient[p,w^2];

"See 'Enlarging quartic ring' script to see rationale behind choices of A, B, F mod a y^2."

m' = m3 - 3/2 r/(2q - b y) m2 + (5q - b y)/6 m1 + 3/4 r (q - b y)/(2q - b y) m0;
p' = CharacteristicPolynomial[m',w];
Coefficient[p',w^4]//Simplify
Coefficient[p',w^3]
Coefficient[p',w^2]//Simplify
"Substitute syzygy in here:"
-(1/(72 (-2 q+b y)^2))y (3 a c (6q(b^2 - 3 a c)y^2 + (9a b c - 2b^3)y^3) y-20 b^4 q y^3+4 b^5 y^4-4 b q (6q(b^2 - 3 a c)y^2 + (9a b c - 2b^3)y^3+18 a c q y^2)+2 b^2 y (6q(b^2 - 3 a c)y^2 + (9a b c - 2b^3)y^3+27 a c q y^2)+12 b^3 (2 q^2 y^2-a c y^4))//Expand
Coefficient[p',w]//Simplify
poly = Coefficient[w p',w]//Simplify;
p3 = Coefficient[poly/19683,r^6]//Simplify
p2 = Coefficient[poly/729,r^4]//Simplify
p1 = Coefficient[poly/27,r^2]//Simplify
p0 = Coefficient[r poly, r]//Simplify

"Now we substitute in syzygy.  Note: Since 3 divides q+by, and we assume 2 divides u, from now on we substitute q by u' = (q+by)/6, and syzygy reads r'^2 = syz"
oldsyz = -8q^3 + 6q(b^2-3a c)y^2 +(9a b c - 2b^3 - 27a^2 d)y^3;
syz = -(u'-b' y)^2 u' + a y^2(-u' c' + (b' c' - a d')y)
coeff0 = p3 oldsyz^3 + p2 oldsyz^2 + p1 oldsyz + p0//Factor;
"Everything is integral!! (Except for some integer denominators...)"

"Now we check that the resolvent map lands in R.  See written notes for details.  Here is just computation.  We need a y^2 to divide A0:"
vec = {W,-u' W - 2a y r' V (1 + n y)/u' + a y^2 U, -r'/(u'-b' y) W + a y V};
A0 = vec[[3]]^2 - vec[[1]] vec[[2]]//Expand//Factor
Coefficient[A0,W^2]//Factor
Coefficient[A0,V W]//Factor
Coefficient[A0,V^2]//Factor
Coefficient[A0,U W]//Factor
Coefficient[A0,U^2]//Factor
Coefficient[A0,U V]//Factor

"Now to work out the omega component.  Need a^2 y^3 to divide B0:"
B0 = ((4b'^2 - 12a c')y^2 - (3u'-2b' y)^2)/3 vec[[1]]^2 + vec[[2]]^2 + (6u'-4b' y) vec[[3]]^2 + 8r' vec[[1]] vec[[3]] - 2u' A0//Expand //Factor
B1 = u'^2 (u' - b' y) B0;
Coefficient[B1,W^2]//Factor
Coefficient[B1,V W]//Factor
Coefficient[B1,V^2]//Factor
Coefficient[B1,U W]//Factor
Coefficient[B1,U^2]//Factor
Coefficient[B1,U V]//Factor

"Everything is integral, yippee! Let's check we actually have a ring now."
basis = {{1,0,3u'/y,8r' (6u'-8b' y)/(4a y^2 (4u'-4b' y))},{0,1,- 2r' (1+n y)/(u' y),(10u'-8b' y)/(2a y^2)},{0,0,1/y,-r'/(a y^2(u'-b' y))},{0,0,0,1/(a y^2)}};
basisinverse = {{1,0,-3u',r'(4b' y-6u')/(u'-b' y)},{0,1,2 r' (1+n y)/u', 2r'^2 (1+n y)/(u'(u'-b' y)) - (5u'-4b' y)},{0,0,y,r' y/(u' - b' y)},{0,0,0,a y^2}};
basis.basisinverse//Simplify
basisinverse.basis//Simplify

"Mult. matrices for powers of alpha (with power basis) with q substituted by u'.  Comparing to non-specialized script, m1new, m2new and m3new are correct:"
m1new = {{0,0,0,((3u'-2b' y)^2-(4b'^2-12a c')y^2)/3},{1,0,0,-8r'},{0,1,0,-6u' + 4b' y},{0,0,1,0}};
m2new = {{0,0,((3u'-2b' y)^2-(4b'^2-12a c')y^2)/3,0},{0,0,-8r',((3u'-2b' y)^2-(4b'^2-12a c')y^2)/3},{1,0,-6u'+4b' y,-8r'},{0,1,0,-6u'+4b' y}};
m3new = {{0,((3u'-2b' y)^2-(4b'^2-12a c')y^2)/3,0,-(6u'-4b' y) ((3u'-2b' y)^2-(4b'^2-12a c')y^2)/3},{0,-8r',((3u'-2b' y)^2-(4b'^2-12a c')y^2)/3,(6u'-4b' y) 8r'},{0,-6u'+4b' y,-8r',(13 (3u'-2b' y)^2-(4b'^2-12a c')y^2)/3},{1,0,-6u'+4b' y,-8r'}};

mbeta = (m2new - 2 r' (1+n y)/u' m1new + 3u' m0)/y;
mgamma = (m3new - r'/(u'-b' y) m2new + (5u'-4b' y) m1new + r' (3u'-4b' y)/(u'-b' y) m0)/(a y^2);

"At a glance, everything looks integral; there are no powers of 2 in denominators."
malpha' = basisinverse.m1new.basis//Expand//Factor
mgamma' = basisinverse.mgamma.basis//Expand//Factor
mbeta' = basisinverse.mbeta.basis//Expand//Factor

(4 a y^4 (Derivative[1][b])^3 Derivative[1][c]-4 y (Derivative[1][b])^2 (-(u'-b' y)^2 u' + a y^2(-u' c' + (b' c' - a d')y))-a y Derivative[1][c] (-(u'-b' y)^2 u' + a y^2(-u' c' + (b' c' - a d')y))-4 y^3 (Derivative[1][b])^4 Derivative[1][u]-9 a y^3 (Derivative[1][b])^2 Derivative[1][c] Derivative[1][u]+4 Derivative[1][b] (-(u'-b' y)^2 u' + a y^2(-u' c' + (b' c' - a d')y)) Derivative[1][u]+12 y^2 (Derivative[1][b])^3 (Derivative[1][u])^2+6 a y^2 Derivative[1][b] Derivative[1][c] (Derivative[1][u])^2-12 y (Derivative[1][b])^2 (Derivative[1][u])^3-a y Derivative[1][c] (Derivative[1][u])^3+4 Derivative[1][b] (Derivative[1][u])^4)//Factor


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
