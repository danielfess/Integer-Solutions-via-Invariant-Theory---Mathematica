#!/usr/bin/env wolframscript
(* ::Package:: *)

"Version with A'', B'', mu = u.  A means A'' and B means B''"

(b^2 y^2/48 - b^3 y^3/144 + u r y^2/4 + u^2 y^4/6 + (3r/2 + 2u y^2) A q y^2 + A^2 y^4 q^2 + B y^2 q^2/3)(b^2 y^2 - 4q^2) + q^5/9 + 2B y^2 q^4/3 + B^2 y ^4 q^3 - (-8q^3 + 6q b^2 y^2 - 2b^3 y^3) q^2/24 - r y^2 q^2 u + B q y^2 (-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/12 + 3B q y^4 r u + 3A B q^2 y^4 r - A r y^2 q^3 //Expand

"Version with original A, B and F substituted by expression in terms of A and B, mod a y^2 - have checked with other Mathematica script that this is correct"
6F^2 + B^2 q - 6A F q + 5A^2 q^2/6 - 7B q^2/3 + 5q^3/4 + 3A B r - 9F r + A q r + (-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/9 + A^2 b^2 y^2/6 + B b^2 y^2/3 - b^2 q y^2/4 //Expand

"Substitute F = A q/2 + 3r/4 mod a y^2"
6(A q/2 + 3r/4)^2 + B^2 q - 6A (A q/2 + 3r/4) q + 5A^2 q^2/6 - 7B q^2/3 + 5q^3/4 + 3A B r - 9(A q/2 + 3r/4) r + A q r + (-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/9 + A^2 b^2 y^2/6 + B b^2 y^2/3 - b^2 q y^2/4 //Expand;

"Result:"
B^2 q-(2 A^2 q^2)/3-(7 B q^2)/3+(13 q^3)/36+3 A B r-(7 A q r)/2-(27 r^2)/8+1/6 A^2 b^2 y^2+1/3 b^2 B y^2+5/12 b^2 q y^2-(2 b^3 y^3)/9

"Substitute r^2 = ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) mod a y^2:"
B^2 q-(2 A^2 q^2)/3-(7 B q^2)/3+(13 q^3)/36+3 A B r-(7 A q r)/2-(-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/8+1/6 A^2 b^2 y^2+1/3 b^2 B y^2+5/12 b^2 q y^2-(2 b^3 y^3)/9//Expand

"Consider this as a quadratic polynomial in A and B.  Leading term is:"
B^2 q-(2 A^2 q^2)/3+3 A B r+1/6 A^2 b^2 y^2

"mod y^2 this factorises as:"
-2/3 (A q - B 9r/(4q))^2

"A factorisation mod a y^2 would have to descend to this so it looks like this, with A ai - B bi = 0 mod y^2:"
-2/3 (A (q+a0) - B (9r/(4q)+b0))(A (q+a1) - B (9r/(4q)+b1))

"The difference between this and the full quadratic term mod a y^2 is:"
-2/3 (A (q+a0) - B (9r/(4q)+b0))(A (q+a1) - B (9r/(4q)+b1)) - (B^2 q-(2 A^2 q^2)/3+3 A B r+1/6 A^2 b^2 y^2) //Expand
-(2/3) A^2 a0 a1+2/3 A a1 B b0+2/3 A a0 B b1-2/3 B^2 b0 b1-2/3 A^2 a0 q-2/3 A^2 a1 q-B^2 q+2/3 A B b0 q+2/3 A B b1 q+(3 A a0 B r)/(2 q)+(3 A a1 B r)/(2 q)-(3 B^2 b0 r)/(2 q)-(3 B^2 b1 r)/(2 q)-(B^2 (-8q^3 + 6q b^2 y^2 - 2b^3 y^3))/(8 q^2)-1/6 A^2 b^2 y^2//Expand

"Integrality of \phi(\gamma) is weaker than integrality of \gamma so let's try to make \phi(\gamma) integral.  This says that the B is congruent to A^2 + q mod a y^2 (theta coeff of phi(gamma)) and some other relation from the omega coeff.
Substituting B == (A^2 + q) mod a y^2 in the full x^2 coeff of min poly we get:"
(A^2 + q)^2 q-(2 A^2 q^2)/3-(7 (A^2 + q) q^2)/3+(49 q^3)/36+3 A (A^2 + q) r-(7 A q r)/2+1/6 A^2 b^2 y^2+1/3 b^2 (A^2 + q) y^2-1/3 b^2 q y^2+(b^3 y^3)/36//Expand

"From quartic form setting, we know that A == -3/4 r q^{-1} solves this mod y (and mod y^2 ??), so we set A = (-3/4 r/q + A' y^2) to find A mod a y^2.  Second step below is reducing 27r^2"
(-3/4 r/q + A' y^2)^4 q-(-3/4 r/q + A' y^2)^2 q^2+q^3/36+3(-3/4 r/q + A' y^2)^3 r-((-3/4 r/q + A' y^2) q r)/2+1/2 (-3/4 r/q + A' y^2)^2 b^2 y^2+(b^3 y^3)/36//Expand
q^3/36-(3 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27))/16-(243 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27)^2)/(256 q^3)+(9 b^2 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y^2)/(32 q^2)+(b^3 y^3)/36+q r y^2 Derivative[1][A]+(27 r ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y^2 Derivative[1][A])/(8 q^2)-(3 b^2 r y^4 Derivative[1][A])/(4 q)-q^2 y^4 (Derivative[1][A])^2-(27 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y^4 (Derivative[1][A])^2)/(8 q)+1/2 b^2 y^6 (Derivative[1][A])^2+q y^8 (Derivative[1][A])^4//Expand

"Make the substitution A' = (br/(q y) A''') and divide by b^4 y^4.  Second calculation is substituting for r^2"
((b^4 y^4)/(64 q)+(b^5 y^5)/(96 q^2)-(b^6 y^6)/(192 q^3)-(b^3 r y^5 (b r/(q y) A'''))/(4 q^2)-1/4 b^2 y^6 (b r/(q y) A''')^2+(b^3 y^7 (b r/(q y) A''')^2)/(4 q)+q y^8 (b r/(q y) A''')^4) q^3/(b^4 y^4) //Expand
q^2/64+(b q y)/96-(b^2 y^2)/192-1/4 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) \!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)-1/4 q ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) (\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\))^2+1/4 b ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y (\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\))^2+((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27)^2 (\!\(\*SuperscriptBox[\(A\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\))^4//Expand

"See notes for justification why A' = br/qy A''' with A''' taking one of the five forms below.  We then substitute each of them in.
Case 1: A''' = (C (q + b y)/(2q - b y)^2). To get rid of denominators I have to compute with homogenised form of polynomial."
(q^2/64+(b q y)/96-(b^2 y^2)/192)(2q - b y)^8-1/4 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) (C (q + b y)(2q - b y)^6)-1/4 q ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) ((C (q + b y)(2q - b y)^2))^2+1/4 b ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y ((C (q + b y)(2q - b y)^2))^2+((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27)^2 ((C (q + b y)))^4//Expand//Factor
"Checking each factor shows no solution for C sadly.
Try case 2: A''' = (C (3q - b y)/(2q - b y)^2).  Same deal with getting rid of denominators."

(q^2/64+(b q y)/96-(b^2 y^2)/192)(2q - b y)^8-1/4 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) (C (3q - b y)(2q - b y)^6)-1/4 q ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) ((C (3q - b y)(2q - b y)^2))^2+1/4 b ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y ((C (3q - b y)(2q - b y)^2))^2+((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27)^2 ((C (3q - b y)))^4//Expand//Factor
"Doesn't work either.
Try case 3: A''' = (C (3q - b y)/((2q - b y)(q + b y))"
(q^2/64+(b q y)/96-(b^2 y^2)/192)(2q - b y)^4 (q + b y)^4 -1/4 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) (C (3q - b y)(2q - b y)^3 (q + b y)^3)-1/4 q ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) ((C (3q - b y)(2q - b y) (q + b y)))^2+1/4 b ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y ((C (3q - b y)(2q - b y) (q + b y)))^2+((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27)^2 ((C (3q - b y)))^4//Expand//Factor
"No luck.  Case 4: A''' = C/(q + b y)"
(q^2/64+(b q y)/96-(b^2 y^2)/192) (q + b y)^4 -1/4 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) (C (q + b y)^3)-1/4 q ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) ((C (q + b y)))^2+1/4 b ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y ((C (q + b y)))^2+((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27)^2 C^4//Expand//Factor
"No luck here either.  Case 5 it is, fingers crossed:  A''' = C/(2q-b y)"
(q^2/64+(b q y)/96-(b^2 y^2)/192) (2q - b y)^4 -1/4 ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) (C (2q - b y)^3)-1/4 q ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) ((C (2q - b y)))^2+1/4 b ((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27) y ((C (2q - b y)))^2+((-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/27)^2 C^4//Expand//Factor

"This works!  C = -3/4 is the unique solution here.  Then using our earlier expression for A (and B in terms of A) and calculating naively i.e. not taking into account whether inverses are taken mod a or mod y^2 or etc..., we get A = (-3/2 r/(2q - b y)), B = ((5q - b y)/6).  Let's test these:"
((5q - b y)/6)^2 q-(2 (-3/2 r/(2q - b y))^2 q^2)/3-(7 ((5q - b y)/6) q^2)/3+(13 q^3)/36+3 (-3/2 r/(2q - b y)) ((5q - b y)/6) r-(7 (-3/2 r/(2q - b y)) q r)/2-(-8q^3 + 6q b^2 y^2 - 2b^3 y^3)/8+1/6 (-3/2 r/(2q - b y))^2 b^2 y^2+1/3 b^2 ((5q - b y)/6) y^2+5/12 b^2 q y^2-(2 b^3 y^3)/9//Expand//Factor

"Yippee!!! This is 0 mod a y^2.  So provided (2q - b y) is invertible mod a y^2 and that we can ignore the denominators in the definitions i.e. we have certain divisibility conditions, then the x^2 coefficient is good to go.  Let's test the rest of the char. poly."




(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
