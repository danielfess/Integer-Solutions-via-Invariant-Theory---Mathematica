#!/usr/bin/env wolframscript
(* ::Package:: *)

a1 = a u^3 - b u^2 s + c u s^2 - d s^3;
b1 = - 3a u^2 t + b (u^2 r + 2 u t s) - c (s^2 t + 2 r s u) + 3d r s^2;
c1 = 3a u t^2 - b (s t^2 + 2 r u t) + c (r^2 u + 2 r s t) - 3d r^2 s;
d1 = - a t^3 + b r t^2 - c r^2 t + d r^3;
l = u r - t s;

c1^2 - 3 b1 d1 - l^2 (r^2 (c^2 - 3b d) + r t(9a d - b c) + t^2(b^2 - 3a c))//Expand;
9 a1 d1 - b1 c1 - l^2 (2r s(c^2 - 3b d) + (r u + t s)(9a d - b c) + 2t u(b^2 - 3a c))//Expand;
b1^2 - 3 a1 c1 - l^2 (s^2(c^2 - 3b d) + s u(9a d - b c) + u^2(b^2 - 3a c))//Expand;

D1 = 9 b c d - 2 c^3 - 27 a d^2;
D2 = 3 b c^2 + 27 a c d - 18 b^2 d;
D3 = 3 b^2 c + 27 a b d - 18 a c^2;
D4 = 9 a b c - 2 b^3 - 27 a^2 d;

9 b1 c1 d1 - 2 c1^3 - 27 a1 d1^2 - l^3 (r^3 D1 + r^2 t D2 + r t^2 D3 + t^3 D4)//Expand;
3 b1 c1^2 + 27 a1 c1 d1 - 18 b1^2 d1 - l^3 (3 r^2 s D1 + (r^2 u + 2 r s t) D2 + (s t^2 + 2 r t u) D3 + 3 t^2 u D4)//Expand;
3 b1^2 c1 + 27 a1 b1 d1 - 18 a1 c1^2 - l^3 (3 r s^2 D1 + (s^2 t + 2 r s u) D2 + (r u^2 + 2 s t u) D3 + 3 t u^2 D4)//Expand;
9 a1 b1 c1 - 2 b1^3 - 27 a1^2 d1 - l^3 (s^3 D1 + s^2 u D2 + s u^2 D3 + u^3 D4)//Expand;
