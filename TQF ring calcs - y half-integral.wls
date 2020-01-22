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
malpha' = basisinverse.m1new.basis//Simplify;
mbeta' = basisinverse.mbeta.basis//Simplify;
mgamma' = basisinverse.mgamma.basis//Simplify;

"Mult matrices above but with syzygy substituted"
malpha'' = {{0,1/2 (-3 u+b y),-((r (b-3 n u+b n y))/(4 u)),-((3 syz +2 (2 u-b y) (6 u^2-3 b u y+a c y^2))/(8 a y^2 (-2 u+b y)))},{1,(r+n r y)/(2 u),(-2 u^2 (-2 u+b y)^2+ syz (1+n y) (b y (1+n y)-u (1+2 n y)))/(4 u^2 y (2 u-b y)),-((r (1+n y) (syz+2 u (-2 u+b y)^2))/(8 a u y^2 (-2 u+b y)^2))},{0,y,r/(4 u-2 b y)-(r+n r y)/(2 u),-((syz+2 u (-2 u+b y)^2)/(4 a y (-2 u+b y)^2))},{0,0,a y,r/(-4 u+2 b y)}}//Simplify;
mbeta'' = {{0,-((r (b-3 n u+b n y))/(4 u)),(2 u^2 (12 u^2-8 b u y+(b^2+a c) y^2)+ syz (1+n y) (u (3-3 n y)+b y (1+n y)))/(8 u^2 y^2),-((r (3 syz (1+n y)+24 u^3 (1+n y)-24 b u^2 y (1+n y)-2 a b c y^3 (1+n y)+2 u y^2 (3 b^2 (1+n y)+a c (3+2 n y))))/(16 a u y^3 (2 u-b y)))},{0,(-2 u^2 (-2 u+b y)^2+ syz (1+n y) (b y (1+n y)-u (1+2 n y)))/(4 u^2 y (2 u-b y)),(r (4 b^2 u^2 y (1+n y)-b (syz (1+n y)^3+8 u^3 (1+2 n y))+2 n u (8 u^3+ syz (1+n y)^2)))/(8 u^3 y (2 u-b y)),(syz^2 (1+n y)^2+4 u^2 (-2 u+b y)^2 (8 u^2-4 b u y+a c y^2)+2 syz u (b^2 y^2 (1+n y)^2-4 b u y (2+2 n y+n^2 y^2)+4 u^2 (3+2 n y+n^2 y^2)))/(16 a u^2 y^3 (-2 u+b y)^2)},{1,r/(4 u-2 b y)-(r+n r y)/(2 u),-((syz (1+n y) (b-2 n u+b n y))/(4 u^2 (2 u-b y))),(r (1+n y) (syz+2 u (-2 u+b y)^2))/(8 a u y^2 (-2 u+b y)^2)},{0,a y,-((a r (1+n y))/u),(2 u (-2 u+b y)^2+ syz (1+n y))/(4 u y (2 u-b y))}}//Simplify;
mgamma'' = {{0,-((3 syz+2 (2 u-b y) (6 u^2-3 b u y+a c y^2))/(8 a y^2 (-2 u+b y))),-((r (3 syz (1+n y)+24 u^3 (1+n y)-24 b u^2 y (1+n y)-2 a b c y^3 (1+n y)+2 u y^2 (3 b^2 (1+n y)+a c (3+2 n y))))/(16 a u y^3 (2 u-b y))),(a c (syz+8 u^3) y+4 b^4 u y^3-4 b u (syz+8 u^3+4 a c u y^2)+2 b^2 y (syz+24 u^3+5 a c u y^2)-2 b^3 (12 u^2 y^2+a c y^4))/(16 a^2 y^3 (-2 u+b y)^2)},{0,-((r (1+n y) (syz+2 u (-2 u+b y)^2))/(8 a u y^2 (-2 u+b y)^2)),(syz^2 (1+n y)^2+4 u^2 (-2 u+b y)^2 (8 u^2-4 b u y+a c y^2)+2 syz u (b^2 y^2 (1+n y)^2-4 b u y (2+2 n y+n^2 y^2)+4 u^2 (3+2 n y+n^2 y^2)))/(16 a u^2 y^3 (-2 u+b y)^2),-((r (b-2 n u+b n y) (syz+(2 u-b y) (4 u^2-2 b u y+a c y^2)))/(8 a^2 u y^3 (-2 u+b y)^2))},{0,-((syz+2 u (-2 u+b y)^2)/(4 a y (-2 u+b y)^2)),(r (1+n y) (syz+2 u (-2 u+b y)^2))/(8 a u y^2 (-2 u+b y)^2),-((syz+(2 u-b y) (4 u^2-2 b u y+a c y^2))/(4 a^2 y^3 (-2 u+b y)))},{1,r/(-4 u+2 b y),(2 u (-2 u+b y)^2+syz (1+n y))/(4 u y (2 u-b y)),0}}//Simplify;

"Now substitute 4b for b, 16c for c, 64d for d, 2y for y1, 4r for r"
malpha''' = {{0,1/2 (-3 u+2b y1),-((4r (4b-3 n u+2b n y1))/(4 u)),(32 c u- 32b c y1+96 a d y1)/(8 (-2 u+2b y1))},{1,(2r+n r y1)/u,(-2 u^2 (-2 u+2b y1)^2+(2+n y1) (- u (-2 u+2b y1)^2-4a y1^2 (c u-b c y1+a d y1)) (b y1 (2+n y1)-u (1+n y1)))/(2 u^2 y1 (2 u-2b y1)),(16r (c u-b c y1+a d y1) (2+n y1))/(2 u (-2 u+2b y1)^2)},{0,y1/2,r/(u- b y1)-(2r+n r y1)/u,(16 y1 (c u-b c y1+a d y1))/(4 (-2 u+2b y1)^2)},{0,0,a y1/2,r/(-u+b y1)}}//Simplify
mbeta''' = {{0,-((r (4b-3 n u+2b n y1))/u),(2 u^2 (12 u^2-16 b u y1+4(b^2+a c) y1^2)+2(2+n y1) (-u (-u+b y1)^2-a y1^2 (c u-b c y1+a d y1)) (u (6-3 n y1)+2b y1 (2+n y1)))/(2 u^2 y1^2),(r (32 c n u-32b (2c+c n y1)+96 a (2d+d n y1)))/(4u (2 u-2b y1))},{0,(-2 u^2 (-2 u+2b y1)^2+(2+n y1) (-u (-2 u+2b y1)^2-4a y1^2 (c u-b c y1+a d y1)) (b y1 (2+n y1)-u (1+n y1)))/(2 u^2 y1 (2 u-2b y1)),(r (16 b^2 u^2 y1 (2+n y1)+2 n u (8 u^3+(2+n y1)^2 (-2 u (-u+b y1)^2-a y1^2 (2 c u-2b c y1+2a d y1)))-4b (8 u^3 (1+n y1)+(2+n y1)^3 (- u (-u+b y1)^2-a y1^2 (c u-b c y1+a d y1)))))/(u^3 y1 (2 u-2b y1)),1/(2 a u^2 y1^3 (-2 u+2b y1)^2) (4 u^2 (-2 u+2b y1)^2 (8 u^2-8 b u y1+4a c y1^2)+(2+n y1)^2 (- u (2 u-2b y1)^2+4a y1^2 (-c u+(b c-a d) y1))^2+2 u (-2 u (2 u-2b y1)^2+8a y1^2 (-c u+(b c-a d) y1)) (b^2 y1^2 (2+n y1)^2-2 b u y1 (8+4n y1+n^2 y1^2)+u^2 (12+4n y1+n^2 y1^2)))},{1,r/(u- b y1)-(2r+n r y1)/u,-(((2+n y1) (2b-n u+b n y1) (-2 u (-2 u+2b y1)^2-8a y1^2 (c u-b c y1+a d y1)))/(4 u^2 (2 u-2b y1))),-((16r (c u-b c y1+a d y1) (2+n y1))/(2u (-2 u+2b y1)^2))},{0,a y1/2,-((a 4r (2+n y1))/(2u)),-((n (u^3-2 b u^2 y1+(b^2+a c) u y1^2+a (-b c+a d) y1^3)+2 a y1 (a d y1+c (u-b y1)))/(u (u-b y1)))}}//Simplify
mgamma''' = {{0,(32 c u-32b c y1+96 a d y1)/(8 (-2 u+2b y1)),(r (32 c n u-32b (2c+c n y1)+96 a (2d+d n y1)))/(4 u (2 u-2b y1)),(-512a c d y1+256c^2 (-2 u+2b y1)-512 b d (-2 u+2b y1))/(16 (-2 u+2b y1)^2)},{0,(16r (c u-b c y1+a d y1) (2+n y1))/(2 u (-2 u+2b y1)^2),1/(2 a u^2 y1^3 (-2 u+2b y1)^2) (4 u^2 (-2 u+2b y1)^2 (8 u^2-8 b u y1+4a c y1^2)+(2+n y1)^2 (4 u (-u+b y1)^2+4a y1^2 (c u-b c y1+a d y1))^2+2 u (-2 u (-2 u+2b y1)^2-8a y1^2 (c u-b c y1+a d y1)) (b^2 y1^2 (2+n y1)^2-2 b u y1 (8+4 n y1+n^2 y1^2)+ u^2 (12+4 n y1+n^2 y1^2))),(64d r (4b-2 n u+2b n y1))/(2 u (-2 u+2b y1)^2)},{0,(16y1 (c u-b c y1+a d y1))/(4 (-2 u+2b y1)^2),-((16r (c u-b c y1+a d y1) (2+n y1))/(2 u (-2 u+2b y1)^2)),(64d)/(-8 u+8 b y1)},{1,r/(-u+ b y1),(2 u (-2 u+2b y1)^2+(2+n y1) (-u (-2 u+2b y1)^2-4a y1^2 (c u-b c y1+a d y1)))/(2 u y1 (2 u-2b y1)),0}}//Simplify


"Check char poly of beta"
poly2 = CharacteristicPolynomial[mbeta',w];
Coefficient[poly2,w^3]//Simplify
Coefficient[poly2,w^2]//Simplify
Coefficient[poly2,w]//Simplify
Coefficient[poly2 w,w]//Simplify


n(a (-b c+a d) y1^3)+2 a y1 (a d y1+c (-b y1))//Factor



(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
