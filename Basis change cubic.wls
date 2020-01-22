#!/usr/bin/env wolframscript
(* ::Package:: *)

g = a x^3 - b x^2 y + c x y^2 - d y^3
tauinvttrans = {{u,-t},{-X/Y',Y/Y'}}
"Above is tau inverse transpose but ignoring the division by det."
gtau = a (u x - t y)^3 - b (u x - t y)^2 (-X/Y' x + Y/Y' y) + c (u x - t y) (-X/Y' x + Y/Y' y)^2 - d (-X/Y' x + Y/Y' y)^3//Expand;
Coefficient[gtau,x^3]//Simplify//Factor
-Coefficient[gtau,x^2 y]//Simplify//Factor
Coefficient[gtau,x y^2]//Simplify//Factor
-Coefficient[gtau,y^3]//Simplify//Factor

"Note: Just read G2 modular forms for this stuff below"
gamma = {{u,-s},{-t,r}}
fgamma = a (u x' - s y')^3 - b (u x' - s y')^2 (-t x' + r y') + c (u x' - s y') (-t x' + r y')^2 - d (-t x' + r y')^3//Expand
a' = Coefficient[fgamma,x'^3]
b' = Coefficient[fgamma,x'^2 y']
c' = Coefficient[fgamma,x' y'^2]
d' = Coefficient[fgamma, y'^3]

"Calculating omega'.theta'"
(u om - t th)(-s om + r th)//Expand
s t a d-r t (-b d - d om - c th)-(-a c - b om - a th) s u+r u a d//Expand
omth = b d r t+ a d s t+a d r u+a c s u+(d r t+b s u)(r om' + t th') +(c r t + a s u)(s om' + u th')//Expand
omconst = Coefficient[omth,om']
thconst = Coefficient[omth,th']
"These tell us how to adjust om' and th' to get a good basis"

"Calculating omega'^2"
(u om - t th)^2//Expand
t^2 (-b d - d om - c th)-2 t u a d + u^2 (-a c - b om - a th)//Expand
om2 = -b d t^2-d (r om' + t th') t^2-2 a d t u-a c u^2-b (r om' + t th') u^2-(a u^2+c t^2)(s om' + u th')//Expand
"Now, we adjust om' so it is part of the good basis.  (om' - const)^2 = om'^2 - 2 const om' + const^2.  Note that coefficients of mult. table match fgamma coefficients."
Coefficient[om2,om']-2thconst - b'//Expand
Coefficient[om2,th'] + a'

"Calculating theta'^2"
(-s om + r th)^2//Expand
(-a c - b om - a th) s^2-2 r s a d+r^2 (-b d - d om - c th)//Expand
th2 = -b d r^2-2 a d r s-a c s^2-(d r^2+b s^2)(r om' + t th')-(c r^2 + a s^2)(s om' + u th')//Expand
"Now, adjust th' as we did with om' to see that coefficients of mult. table match fgamma coefficients:"
Coefficient[th2,th'] - 2omconst + c'//Expand
Coefficient[th2,om'] - d'






