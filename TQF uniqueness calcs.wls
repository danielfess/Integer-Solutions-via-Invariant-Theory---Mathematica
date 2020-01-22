#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]
"Coefficients of phi(gamma_3) having removed Y^2 factors"
p6 = -3A^2(b^2+(a b c-a^2 d)Y v + a^2 b d Y^2 v^2) - 5A B b W v + 2B^2 (Z+b Y);
q = -A^2(b+a c Y v + a^2 d Y^2 v^2) - 2A B W v + B^2 Y;
r = -a A^2;

"Coefficients of char poly of phi(gamma_3), with factor of Y^2r omitted from associated coeff"
"Double x^2 coeff"
x2coeff = -p6 + 2b q + 2c r//Simplify;
"12 times x coeff"
xcoeff = p6^2 - 4c p6 r - 4b p6 q + 12(b c - 3a d)q r + 12a c q^2 + 12b d r^2//Simplify;
Coefficient[xcoeff, A^4]//Expand;
Coefficient[xcoeff, A^3 B]//Expand;
Coefficient[xcoeff, A^2 B^2]//Expand;
Coefficient[xcoeff, A B^3]//Expand;
Coefficient[xcoeff, B^4]//Expand;

"Not sure what calculations above are.  Repurposing this script for calculations to go with LaTeX document 'Towards TQF Syzygy Uniqueness'."
"Calculating \phi(\lambda_1 \lambda_2 \lambda_3 \beta_3) below"
Bf = {{s,0,4W'},{0,1,0},{4W',0,2Z'}};
A0 = {{0,-1/2,0},{-1/2,0,0},{0,0,1}};
"Vector of lambda1 lambda2 lambda3 beta3 in correct basis of Z[\alpha] for computing resolvent map:"
"G in place of E because for some reason Mathematica writes 'E' as 'e' i.e. exponential?"
lambeta3 = {{l1 l2, B D + G l2 - 2Z' l1 l2, D l1}};
"li = lambda_i"
philambdabeta3om = lambeta3.Bf.Transpose[lambeta3] Y + lambeta3.A0.Transpose[lambeta3] ((4Y Z'-b Y^2)/3 - 2Z' Y)//Expand
philambdabeta3th = -lambeta3.A0.Transpose[lambeta3] a Y^2//Expand

"Calculating multiplication matrices in basis of Q"
"First, via basis of Z[alpha] given by powers of alpha"

m0 = IdentityMatrix[4];
m1 = {{0,0,0,(4Z'^2-(b^2-3a c)y^2)/12},{1,0,0,-8W'},{0,1,0,-2Z'},{0,0,1,0}};
m2 = {{0,0,(4Z'^2-(b^2-3a c)y^2)/12,0},{0,0,-8W',(4Z'^2-(b^2-3a c)y^2)/12},{1,0,-2Z',-8W'},{0,1,0,-2Z'}};
m3 = {{0,(4Z'^2-(b^2-3a c)y^2)/12,0,-2Z' (4Z'^2-(b^2-3a c)y^2)/12},{0,-8W',(4Z'^2-(b^2-3a c)y^2)/12,16 Z' W'},{0,-2Z',-8W',(52 Z'^2-(b^2-3a c)y^2)/12},{1,0,-2Z',-8W'}};

T = Transpose[{{1,0,0,0},{A/l1,1/l1,0,0},{(C+A B')/l2,B'/l2,1/l2,0},{F/l3+A G/(l1 l3) + (A B'+C) D/(l2 l3),G/(l1 l3) + B' D/(l2 l3),D/(l2 l3),1/l3}}]
"Matrices o beta_i in new basis with denominators"
mb1 = Inverse[T].(m1+A m0).T/l1 //Expand;
mb2 = Inverse[T].(m2+C m0).T/l2 + l1 B' mb1/l2//Expand;
mb3 = Inverse[T].(m3+F m0).T/l3 + D mb2/l3 + G mb1/l3 //Expand;
"Mult matrices of beta_i in new basis, rescaled as follows:
beta1 times l1^2 l2 l3
beta2 times l1 l2^2 l3
beta3 times l1^2 l2^2 l3^2"
mb1new = l1^2 l2 l3 mb1//Expand;
mb2new = l1 l2^2 l3 mb2//Expand;
mb3new = l1^2 l2^2 l3^2 mb3//Expand;

"Replace D by l2 D2, G by l1 G':  [Note: Functions renamed with a D because I spotted divisibility with D before that of G.]"
TD = Transpose[{{1,0,0,0},{A/l1,1/l1,0,0},{(C+A B')/l2,B'/l2,1/l2,0},{F/l3+A G'/l3 + (A B'+C)D2/l3,G'/l3 + B' D2/l3,D2/l3,1/l3}}]
"Matrices o beta_i in new basis with denominators"
mb1D = Inverse[TD].(m1+A m0).TD/l1 //Expand;
mb2D = Inverse[TD].(m2+C m0).TD/l2 + l1 B' mb1D/l2//Expand;
mb3D = Inverse[TD].(m3+F m0).TD/l3 + l2 D2 mb2D/l3 + l1 G' mb1D/l3 //Expand;
"Mult matrices of beta_i in new basis, rescaled as follows:
beta1 times l1^2 l2 l3
beta2 times l1 l2^2 l3
beta3 times l1^2 l2^2 l3^2"
mb1Dnew = l1^2 l2 l3 mb1D//Expand
mb2Dnew = l1 l2^2 l3 mb2D//Expand
mb3Dnew = l1^2 l2^2 l3^2 mb3D//Expand






