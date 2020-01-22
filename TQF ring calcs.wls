#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]

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

"Now we substitute in syzygy.  Note: Since 3 divides q+by, from now on we substitute q by u = (q+by)/3, and syzygy reads r^2 = syz"
oldsyz = -8q^3 + 6q(b^2-3a c)y^2 +(9a b c - 2b^3 - 27a^2 d)y^3;
syz = -2(2u-b y)^2 u + a y^2(-2u c + (b c - a d)y);
coeff0 = p3 oldsyz^3 + p2 oldsyz^2 + p1 oldsyz + p0//Factor
"Everything is integral!! (Except for some integer denominators...)"

"Now we check that the resolvent map lands in R.  See written notes for details.  Here is just computation.  We need a y^2 to divide A0:"
vec = {W,-u W/2 - 1/2 a y r V (1 + n y)/u + a y^2 U, -1/2 r/(2u-b y) W + a y V}
A0 = vec[[3]]^2 - vec[[1]] vec[[2]]//Expand//Factor
A1 = 4u (2u-b y)^2 A0;
Coefficient[A1,W^2]//Factor
Coefficient[A1,V W]//Factor
Coefficient[A1,V^2]//Factor
Coefficient[A1,U W]//Factor
Coefficient[A1,U^2]//Factor
Coefficient[A1,U V]//Factor

"Now to work out the omega component.  Need a^2 y^3 to divide B0:"
B0 = ((b^2 - 3a c)y^2 - (3u-b y)^2)/12 vec[[1]]^2 + vec[[2]]^2 + (3u-b y) vec[[3]]^2 + r vec[[1]] vec[[3]] - u A0//Expand //Factor
B1 = 4u^2 (2u - b y) B0;
Coefficient[B1,W^2]//Factor
Coefficient[B1,V W]//Factor
Coefficient[B1,V^2]//Factor
Coefficient[B1,U W]//Factor
Coefficient[B1,U^2]//Factor
Coefficient[B1,U V]//Factor

"Missing some factors above... let's check we actually have a ring first.  Edit: alpha term of beta was wrong, I think.  Should be fixed now."
basis = {{1,0,(3u-b y)/(2y),r (3u-2b y)/(4a y^2 (2u-b y))},{0,1,-1/2 r (1+n y)/(u y),(5u-2b y)/(2a y^2)},{0,0,1/y,-r/(2a y^2(2u-b y))},{0,0,0,1/(a y^2)}};
basisinverse = {{1,0,-(3u-b y)/2,-3r/4},{0,1,1/2 r (1+n y)/u, r^2 (1+n y)/(4u(2u-b y)) - (5u-2b y)/2},{0,0,y,r y/(4u - 2b y)},{0,0,0,a y^2}};
basis.basisinverse//Simplify;
basisinverse.basis//Simplify;

"Mult. matrices for powers of alpha (with basis = power basis) with q substituted by e:"
m1new = {{0,0,0,((3u-b y)^2-(b^2-3a c)y^2)/12},{1,0,0,-r},{0,1,0,-3u + b y},{0,0,1,0}};
m2new = {{0,0,((3u-b y)^2-(b^2-3a c)y^2)/12,0},{0,0,-r,((3u-b y)^2-(b^2-3a c)y^2)/12},{1,0,-3u+b y,-r},{0,1,0,-3u+b y}};
m3new = {{0,((3u-b y)^2-(b^2-3a c)y^2)/12,0,-(3u-b y) ((3u-b y)^2-(b^2-3a c)y^2)/12},{0,-r,((3u-b y)^2-(b^2-3a c)y^2)/12,(3u-b y) r},{0,-3u+b y,-r,(13 (3u-b y)^2-(b^2-3a c)y^2)/12},{1,0,-3u+b y,-r}};

mbeta = (m2new - 1/2 r (1+n y)/u m1new + (3u-b y)/2 m0)/y;
mgamma = (m3new - r/(2(2u-b y)) m2new + (5u-2b y)/2 m1new + r (3u-2b y)/(4(2u-b y)) m0)/(a y^2);

"Note that the only obstruction to integrality below is by powers of 4"
malpha' = basisinverse.m1new.basis//Simplify
mbeta' = basisinverse.mbeta.basis//Simplify
mgamma' = basisinverse.mgamma.basis//Simplify

poly2 = CharacteristicPolynomial[mbeta',w];
Coefficient[poly2,w^3]//Simplify
Coefficient[poly2,w^2]//Simplify
Coefficient[poly2,w]//Simplify
Coefficient[poly2 w,w]//Simplify

8 u^3+(-8u^3+8u^2 b y-2u b^2 y^2+a y^2(-2u c+(b c-a d)y))(1+n y)^2-4 b u^2 y//Factor
-4 b u^2+16 n u^3-16 b n u^2 y+8 n^2 u^3 y-8 b n^2 u^2 y^2
2 b^2 u y+2 a c u y-a b c y^2+a^2 d y^2+4 b^2 n u y^2+4 a c n u y^2-2 a b c n y^3+2 a^2 d n y^3+2 b^2 n^2 u y^3+2 a c n^2 u y^3-a b c n^2 y^4+a^2 d n^2 y^4//Factor

"Checking resulting TQFs solve the system. Because I have used A, B already, omegaform is B and thetaform is A.  First check invariants. This is new result for general I."
thetaform = {{0,0,1/2},{0,-a,0},{1/2,0,-(W^2 g^2+8Z)/(a y^2)}}
omegaform = {{y/16,-W g/8,0},{-W g/8,(W^2 g^2+8Z-4b y)/(4y),W (2-(2Z-b y)g)/(a y^2)},{0,W (2-(2Z-b y)g)/(a y^2),((W^2 g^2-8Z)(2Z-b y)-4W^2 g-4a c y^2)/(a^2 y^3)}}

p = 4 Det[thetaform X + omegaform Y]
Coefficient[p,X^3]
Coefficient[p,X^2 Y]//Simplify
Coefficient[p,X Y^2]//Simplify
Coefficient[p,Y^3]//Simplify

"Checking resulting TQFs solve the system. This is original result."
thetaform2 = {{0,0,1/2},{0,-a,0},{1/2,0,-(W'^2 f^2+Z')/(a y^2)}}
omegaform2 = {{y,-2W' f,0},{-2W' f,4(W'^2 f^2+Z'-b y)/y,4W'(1-(Z'-b y)f)/(a y^2)},{0,4W'(1-(Z'-b y)f)/(a y^2),4((W'^2 f^2-Z')(Z'-b y)-2W'^2 f-a c y^2)/(a^2 y^3)}}

qform = 4 Det[thetaform2 X + omegaform2 Y]
Coefficient[qform,X^3]
Coefficient[qform,X^2 Y]//Simplify
Coefficient[qform,X Y^2]//Simplify
Coefficient[qform,Y^3]//Simplify

"From 'Ternary QFs script', -8j = r, with j[-1,0,0] given in this script leading to simple formula for r[-1,0,0], using a_11=a_12=a_23=b_13=0:"
rfunction= (1/8) (\!\(
\*SubsuperscriptBox[\(a\), \(13\), \(3\)]\ 
\*SubsuperscriptBox[\(b\), \(12\), \(3\)]\)-4 \!\(\*SubsuperscriptBox[\(a\), \(13\), \(3\)]\) Subscript[b, 11] Subscript[b, 12] Subscript[b, 22]-8 \!\(\*SubsuperscriptBox[\(a\), \(13\), \(2\)]\) Subscript[a, 22] \!\(\*SubsuperscriptBox[\(b\), \(11\), \(2\)]\) Subscript[b, 23])
"Evaluating for original result:"
(1/8) ((-4W' f)^3-16 y (-4W' f)(W'^2 f^2+Z'-b y)/y +64a y^2 W'(1-(Z'-b y)f)/(a y^2))
"= 8W'"

"Evaluating for new result.  Note: Need to multiply by 64 because r[-4,0,0] = 64 r[-1,0,0]"
8((-W g/4)^3-4(y/16)(-W g/4)(W^2 g^2+8Z-4b y)/(4y)  -8 (-a)(y/16)^2 (2W (-g(2Z-b y)+2)/(a y^2)))
"Result is W"

"Checking q covariant is correct.  q = -1/2 tr(gamma^2)=-1/2 tp[x,y,z,2]
qfunction is q evaluated at (-1,0,0) using a_11=a_12=a_23=b_13=0"
qfunction = -1/2 (-4 Subscript[a, 22] Subscript[a, 33] \!\(\*SubsuperscriptBox[\(b\), \(11\), \(2\)]\)+3/4 \!\(\*SubsuperscriptBox[\(a\), \(13\), \(2\)]\) \!\(\*SubsuperscriptBox[\(b\), \(12\), \(2\)]\)-2 \!\(\*SubsuperscriptBox[\(a\), \(13\), \(2\)]\) Subscript[b, 11] Subscript[b, 22])
"For original result, we get 6Z'-4b y = q:"
-1/2 (-4 a (W'^2 f^2+Z')/(a y^2) y^2+3/4 (-4W' f)^2-8 y (W'^2 f^2+Z'-b y)/y)//Simplify
"For new result, q(-4,0,0) = 16 q(-1,0,0) = 3Z-by is correct"
-8 (-4 a (W^2 g^2+8Z)/(a y^2) (y/16)^2+3/4 (-W g/4)^2-2 y/16 (W^2 g^2+8Z-4b y)/(4y))//Simplify
"So all covariants are correct!! Yippee!!"

"Work on unicity of solution.  diff = a Y^3(gamma1 - gamma2), where gamma1/2 are the final elements in a basis of Q/Q' (resp.) and we want to see that Q = Q'."
mdiff = C1 m2 + (-3r t zed C1 + C2 y^3) m1 + t q C1 m0;
"where t = inverse of 2 mod a Y^3, and zed = inverse of Z mod Y^3"
"Considering a.diff we need C1/Y^2 beta to be integral, so let's consider the char poly of Y.beta (with q instead of u)...
...actually, in the first place, computation of char poly of beta goes wrong so let's derive from first principles what beta is:"
mbetaq = m2 + A m1 + t q m0;
charbeta = CharacteristicPolynomial[mbetaq,w];
Coefficient[charbeta,w^4];
Coefficient[charbeta,w^3]//Simplify;
Coefficient[charbeta,w^2]//Simplify;
"From these we can see that the only possible choice (which may or may not work) is A = -3W zed t.
We now see for which A is A/Y^2 beta = A/Y^3 (alpha^2+...) integral, as this was our original question:"
mbetaq2 = m2 -3 r zed t m1 + t q m0;
charbeta2 = CharacteristicPolynomial[mbetaq2,w];
Coefficient[charbeta2,w^4]
Coefficient[charbeta2,w^3]//Simplify
Coefficient[charbeta2,w^2]//Simplify
Coefficient[charbeta2,w]//Simplify
Coefficient[w charbeta2,w]//Simplify

chardiff = CharacteristicPolynomial[mdiff,w];
Coefficient[chardiff,w^4]
Coefficient[chardiff,w^3]//Simplify
Coefficient[chardiff,w^2]//Simplify
Coefficient[chardiff,w]//Simplify
Coefficient[w chardiff,w]//Simplify



(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
