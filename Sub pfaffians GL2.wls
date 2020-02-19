#!/usr/bin/env wolframscript
(* ::Package:: *)

g = {{a^3, a^2 b, a b^2, b^3},{3a^2 c, 2a b c + a^2 d, 2a b d + b^2 c, 3b^2 d},{3a c^2, 2a c d + b c^2, 2b c d + a d^2, 3b d^2},{c^3, c^2 d, c d^2, d^3}}

Q1 = {{0,0,1/2,0},{0,-1,0,0},{1/2,0,0,0},{0,0,0,0}}
Q2 = {{0,0,0,0},{0,0,0,1/2},{0,0,-1,0},{0,1/2,0,0}}
Q3 = {{0,0,0,1/2},{0,0,-1/2,0},{0,-1/2,0,0},{1/2,0,0,0}}
Q4f = {{0,f0/2,0,0},{f0/2,f1,f2/2,0},{0,f2/2,f3,f4/2},{0,0,f4/2,f5}}
Q5f = {{f0,f1/2,0,0},{f1/2,f2,f3/2,0},{0,f3/2,f4,f5/2},{0,0,f5/2,0}}

x1 = a x + c y
y1 = b x + d y

G = f0 x1^5 + f1 x1^4 y1 + f2 x1^3 y1^2 + f3 x1^2 y1^3 + f4 x1 y1^4 + f5 y1^5
g0 = Coefficient[G,x^5]
g1 = Coefficient[G,x^4 y]//Expand
g2 = Coefficient[G,x^3 y^2]//Expand
g3 = Coefficient[G,x^2 y^3]//Expand
g4 = Coefficient[G,x y^4]//Expand
g5 = Coefficient[G,y^5]

Q4g = {{0,g0/2,0,0},{g0/2,g1,g2/2,0},{0,g2/2,g3,g4/2},{0,0,g4/2,g5}}
Q5g = {{g0,g1/2,0,0},{g1/2,g2,g3/2,0},{0,g3/2,g4,g5/2},{0,0,g5/2,0}}

g.Q1.Transpose[g] - (a d - b c)^2 (a^2 Q1 + c^2 Q2 + a c Q3)//Simplify
g.Q2.Transpose[g] - (a d - b c)^2 (b^2 Q1 + d^2 Q2 + b d Q3)//Simplify
g.Q3.Transpose[g] - (a d - b c)^2 (2a b Q1 + 2c d Q2 + (b c + a d) Q3)//Simplify
g.Q4f.Transpose[g] - b Q5g - d Q4g - (3 a^3 b c^2 f0+a c (b c+2 a d) (a^2 f0+2 a b f1+b^2 f2)+b d (2 b c+a d) (a^2 f2+2 a b f3+b^2 f4)+3 b^3 d^2 (a f4+2 b f5))Q1 - (a^3 c^2 d f0+a^2 b c (c^2 f0+2 c d f1+d^2 f2)+a b^2 d (c^2 f2+2 c d f3+d^2 f4)+b^3 d^2 (c f4+2 d f5))Q3 - (a^2 c d (4 c^2 f0+2 c d f1+d^2 f2)+2 a b (c^4 f0+2 c^3 d f1+2 c^2 d^2 f2+2 c d^3 f3+d^4 f4)+b^2 d (c^3 f2+2 c^2 d f3+4 c d^2 f4+6 d^3 f5))Q2//Simplify

(3 a^3 b c^2 f0+a c (b c+2 a d) (a^2 f0+2 a b f1+b^2 f2)+b d (2 b c+a d) (a^2 f2+2 a b f3+b^2 f4)+3 b^3 d^2 (a f4+2 b f5))//Simplify
(a^3 c^2 d f0+a^2 b c (c^2 f0+2 c d f1+d^2 f2)+a b^2 d (c^2 f2+2 c d f3+d^2 f4)+b^3 d^2 (c f4+2 d f5))//Simplify
(a^2 c d (4 c^2 f0+2 c d f1+d^2 f2)+2 a b (c^4 f0+2 c^3 d f1+2 c^2 d^2 f2+2 c d^3 f3+d^4 f4)+b^2 d (c^3 f2+2 c^2 d f3+4 c d^2 f4+6 d^3 f5))//Simplify

g.Q5f.Transpose[g] - a Q5g - c Q4g - (3 a^3 c^2 (2 a f0+b f1)+a c (b c+2 a d) (a^2 f1+2 a b f2+b^2 f3)+3 a b^3 d^2 f5+b d (2 b c+a d) (a^2 f3+2 a b f4+b^2 f5))Q1 - (a^2 c (6 c^3 f0+4 c^2 d f1+2 c d^2 f2+d^3 f3)+b^2 c d (c^2 f3+2 c d f4+4 d^2 f5)+2 a b (c^4 f1+2 c^3 d f2+2 c^2 d^2 f3+2 c d^3 f4+d^4 f5))Q2 - (a^3 c^2 (2 c f0+d f1)+a^2 b c (c^2 f1+2 c d f2+d^2 f3)+b^3 c d^2 f5+a b^2 d (c^2 f3+2 c d f4+d^2 f5))Q3//Simplify
(3 a^3 c^2 (2 a f0+b f1)+a c (b c+2 a d) (a^2 f1+2 a b f2+b^2 f3)+3 a b^3 d^2 f5+b d (2 b c+a d) (a^2 f3+2 a b f4+b^2 f5))//Simplify
(a^2 c (6 c^3 f0+4 c^2 d f1+2 c d^2 f2+d^3 f3)+b^2 c d (c^2 f3+2 c d f4+4 d^2 f5)+2 a b (c^4 f1+2 c^3 d f2+2 c^2 d^2 f3+2 c d^3 f4+d^4 f5))//Simplify
(a^3 c^2 (2 c f0+d f1)+a^2 b c (c^2 f1+2 c d f2+d^2 f3)+b^3 c d^2 f5+a b^2 d (c^2 f3+2 c d f4+d^2 f5))//Simplify

A = {{0, -f2 x - f3 y, -f0 w - f1 x, y, z},{f2 x + f3 y, 0, f4 y + f5 z, w, x},{f0 w + f1 x, - f4 y - f5 z, 0, -x, -y},{- y, - w, x, 0, 0},{- z, - x, y, 0, 0}}
subpfaff = Reverse[Diagonal[Minors[A]]]//Factor




