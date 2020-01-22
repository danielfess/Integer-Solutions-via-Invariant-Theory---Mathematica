#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]
m0 = IdentityMatrix[4]
m1 = {{0,3x^2 - 4i y^2,0,-8z},{0,0,0,1},{1,0,0,6x},{0,0,1,0}}
m2 = {{0,0,-8z,3x^2 - 4i y^2},{0,0,1,0},{0,3x^2 - 4i y^2,6x,-8z},{1,0,0,6x}}
m3 = {{-8z,0,3x^2 - 4i y^2,-48x z},{1,0,0,6x},{6x,0,-8z,39x^2 - 4i y^2},{0,3x^2 - 4i y^2,6x,-8z}}
c = (1 + 2e) z/x;
m = m3 + (c+2 z/x)m2 + (-6x + c^2 + y f2)m1 + (-3x c + y f2 c)m0;
m//MatrixForm;
p = CharacteristicPolynomial[m,w]//Expand;
Coefficient[p,w^4];
Coefficient[p,w^3]
Coefficient[p,w^2]

"From here, looking at w^3 coeff gives y divides f2 and subsequently from w^2 coeff we see e=-1"

c' = -z/x;
m' = m3 + z/x m2 + (-5x - (4i x y^2 + 4j y^3)/(4 x^2) + y^2 g2)m1 + (3z - y^2 g2 z/x)m0;
p' = CharacteristicPolynomial[m',w]//Expand;
Coefficient[p',w^4];
Coefficient[p',w^3];
Coefficient[p',w^2];
Coefficient[p',w];
Coefficient[w p',w];

basis = {{1,-5z/y^2,0,-3x/y},{0,1/y^2,0,0},{0,x/y^2,1,2 z/(x y)},{0,z/(x y^2),0,1/y}};
basisinverse = {{1,5z-3z,0,3x},{0,y^2,0,0},{0,x- (2i x y^2 + 2j y^3)/x^2,1,-2z/x},{0,-z y/x,0,y}};
malpha1 = 1/y^2 (m3 + z/x m2 - 5x m1 + 3z m0);
malpha3 = 1/y (m2 + 2z f m1 - 3x m0);
"f is x^{-1}"
pbeta = CharacteristicPolynomial[malpha3,w];
Coefficient[pbeta,w^4]
Coefficient[pbeta,w^3]
Coefficient[pbeta,w^2]
Coefficient[pbeta,w]
Coefficient[pbeta w,w]

malpha2 = m1;

malpha1' = basisinverse.malpha1.basis//Expand//MatrixForm
malpha2' = basisinverse.malpha2.basis//Expand//MatrixForm
malpha3' = basisinverse.malpha3.basis//Expand//MatrixForm





(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
