#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]
g = f0 (a x + c y)^4 + f1 (a x + c y)^3 (b x + d y) + f2 (a x + c y)^2 (b x + d y)^2 + f3 (a x + c y) (b x + d y)^3 + f4 (b x + d y)^4//Expand
g0 = Coefficient[g,x^4]
g1 = Coefficient[g,x^3 y]
g2 = Coefficient[g,x^2 y^2]
g3 = Coefficient[g,x y^3]
g4 = Coefficient[g,y^4]

P= {{a^2, 2a c, c^2},{a b, a d + b c, c d},{b^2, 2b d, d^2}}
Bf = {{f0,f1/2,0},{f1/2,f2,f3/2},{0,f3/2,f4}}
Bg = {{g0,g1/2,0},{g1/2,g2,g3/2},{0,g3/2,g4}}
Transpose[P].Bf.P - Bg//Expand



