#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]

"Conjugates of phi ignoring factor of 1/sqrt(disc)"

phi1 = x1 x2 + x2 x3 + x3 x4 + x4 x5 + x5 x1 - x1 x3 - x3 x5 - x5 x2 - x2 x4 - x4 x1;
phi2 = x3 x5 + x5 x2 + x2 x1 + x1 x4 + x4 x3 - x3 x2 - x2 x4 - x4 x5 - x5 x1 - x1 x3;
phi3 = x4 x1 + x1 x3 + x3 x2 + x2 x5 + x5 x4 - x4 x3 - x3 x5 - x5 x1 - x1 x2 - x2 x4;
phi4 = x5 x2 + x2 x4 + x4 x3 + x3 x1 + x1 x5 - x5 x4 - x4 x1 - x1 x2 - x2 x3 - x3 x5;
phi5 = x1 x3 + x3 x5 + x5 x4 + x4 x2 + x2 x1 - x1 x5 - x5 x2 - x2 x3 - x3 x4 - x4 x1;
phi6 = x2 x4 + x4 x1 + x1 x5 + x5 x3 + x3 x2 - x2 x1 - x1 x3 - x3 x4 - x4 x5 - x5 x2;

p1 = SymmetricPolynomial[1,{phi1,phi2,phi3,phi4,phi5,phi6}];
p2 = SymmetricReduction[SymmetricPolynomial[2,{phi1,phi2,phi3,phi4,phi5,phi6}],{x1,x2,x3,x4,x5},{s1,s2,s3,s4,s5}];
p3 = SymmetricReduction[SymmetricPolynomial[3,{phi1,phi2,phi3,phi4,phi5,phi6}],{x1,x2,x3,x4,x5},{s1,s2,s3,s4,s5}];
p4 = SymmetricReduction[SymmetricPolynomial[4,{phi1,phi2,phi3,phi4,phi5,phi6}],{x1,x2,x3,x4,x5},{s1,s2,s3,s4,s5}];
p5 = SymmetricReduction[SymmetricPolynomial[5,{phi1,phi2,phi3,phi4,phi5,phi6}],{x1,x2,x3,x4,x5},{s1,s2,s3,s4,s5}];
p6 = SymmetricReduction[SymmetricPolynomial[6,{phi1,phi2,phi3,phi4,phi5,phi6}],{x1,x2,x3,x4,x5},{s1,s2,s3,s4,s5}];
"Note p1,p3 = 0"

"Fixing 5th sym poly because of factor of sqrt(disc)"
sqdiscx = (x1-x2) (x1-x3) (x2-x3) (x1-x4) (x2-x4) (x3-x4) (x1-x5) (x2-x5) (x3-x5) (x4-x5);
p5' = SymmetricReduction[SymmetricPolynomial[5,{phi1,phi2,phi3,phi4,phi5,phi6}] sqdiscx,{x1,x2,x3,x4,x5},{s1,s2,s3,s4,s5}];
discx = Discriminant[w^5-s1 w^4 + s2 w^3 - s3 w^2 + s4 w - s5,w];
p5'[[1]]-32discx//Expand;
"so p5'[[1]] = 32d"

"Adding in 1/disc factors, small tweak needed for 5th coeff as handled above.
However, this is slow computationally, so we add the required factor in afterwards."
qgen = w^6 + S2 w^4 + S4 w^2 - S5 w^1 + S6;
discqgen = Discriminant[qgen,w];
16 S2^4 S4^3 S5^2-128 S2^2 S4^4 S5^2+256 S4^5 S5^2+108 S2^5 S5^4-900 S2^3 S4 S5^4+2000 S2 S4^2 S5^4+3125 S5^6-64 S2^4 S4^4 S6+512 S2^2 S4^5 S6-1024 S4^6 S6-576 S2^5 S4 S5^2 S6+4816 S2^3 S4^2 S5^2 S6-10560 S2 S4^3 S5^2 S6+1500 S2^2 S5^4 S6-22500 S4 S5^4 S6+512 S2^5 S4^2 S6^2-4352 S2^3 S4^3 S6^2+9216 S2 S4^4 S6^2-192 S2^4 S5^2 S6^2-6480 S2^2 S4 S5^2 S6^2+43200 S4^2 S5^2 S6^2-1024 S2^6 S6^3+9216 S2^4 S4 S6^3-17280 S2^2 S4^2 S6^3-13824 S4^3 S6^3-32400 S2 S5^2 S6^3-13824 S2^3 S6^4+62208 S2 S4 S6^4-46656 S6^5
"Note that whenever S5 appears, it appears to an even power.  This is to preserve homoegeneity of disc, since the other odd coeffs are 0."

q2 = p2[[1]]/discx;
q4 = p4[[1]]/discx^2;
q5 = 32/discx^2;
q5sq = 32^2/discx^4;
q6 = p6[[1]]/discx^3;
q = w^6 + q2 w^4 + q4 w^2 - q5 w + q6;
"Want to calculate discq = discphi, but it is computationally slow, so we calculate with the below polys:"
Q2 = p2[[1]];
Q4 = p4[[1]];
Q5sq = 32^2 discx;
Q6 = p6[[1]];
"Q = w^6 + Q2 w^4 + Q4 w^2 - Q5 w^1 + Q6 has roots which are sqdiscx times the roots of q, hence disc(Q) = discx^15 . disc(q)"
"discQ =  16 \!\(\*SuperscriptBox[\(Q2\), \(4\)]\) \!\(\*SuperscriptBox[\(Q4\), \(3\)]\) Q5sq-128 \!\(\*SuperscriptBox[\(Q2\), \(2\)]\) \!\(\*SuperscriptBox[\(Q4\), \(4\)]\) Q5sq+256 \!\(\*SuperscriptBox[\(Q4\), \(5\)]\) Q5sq+108 \!\(\*SuperscriptBox[\(Q2\), \(5\)]\) \!\(\*SuperscriptBox[\(Q5sq\), \(2\)]\)-900 \!\(\*SuperscriptBox[\(Q2\), \(3\)]\) Q4 \!\(\*SuperscriptBox[\(Q5sq\), \(2\)]\)+2000 Q2 \!\(\*SuperscriptBox[\(Q4\), \(2\)]\) \!\(\*SuperscriptBox[\(Q5sq\), \(2\)]\)+3125 \!\(\*SuperscriptBox[\(Q5sq\), \(3\)]\)-64 \!\(\*SuperscriptBox[\(Q2\), \(4\)]\) \!\(\*SuperscriptBox[\(Q4\), \(4\)]\) Q6+512 \!\(\*SuperscriptBox[\(Q2\), \(2\)]\) \!\(\*SuperscriptBox[\(Q4\), \(5\)]\) Q6-1024 \!\(\*SuperscriptBox[\(Q4\), \(6\)]\) Q6-576 \!\(\*SuperscriptBox[\(Q2\), \(5\)]\) Q4 Q5sq Q6+4816 \!\(\*SuperscriptBox[\(Q2\), \(3\)]\) \!\(\*SuperscriptBox[\(Q4\), \(2\)]\) Q5sq Q6-10560 Q2 \!\(\*SuperscriptBox[\(Q4\), \(3\)]\) Q5sq Q6+1500 \!\(\*SuperscriptBox[\(Q2\), \(2\)]\) \!\(\*SuperscriptBox[\(Q5sq\), \(2\)]\) Q6-22500 Q4 \!\(\*SuperscriptBox[\(Q5sq\), \(2\)]\) Q6+512 \!\(\*SuperscriptBox[\(Q2\), \(5\)]\) \!\(\*SuperscriptBox[\(Q4\), \(2\)]\) \!\(\*SuperscriptBox[\(Q6\), \(2\)]\)-4352 \!\(\*SuperscriptBox[\(Q2\), \(3\)]\) \!\(\*SuperscriptBox[\(Q4\), \(3\)]\) \!\(\*SuperscriptBox[\(Q6\), \(2\)]\)+9216 Q2 \!\(\*SuperscriptBox[\(Q4\), \(4\)]\) \!\(\*SuperscriptBox[\(Q6\), \(2\)]\)-192 \!\(\*SuperscriptBox[\(Q2\), \(4\)]\) Q5sq \!\(\*SuperscriptBox[\(Q6\), \(2\)]\)-6480 \!\(\*SuperscriptBox[\(Q2\), \(2\)]\) Q4 Q5sq \!\(\*SuperscriptBox[\(Q6\), \(2\)]\)+43200 \!\(\*SuperscriptBox[\(Q4\), \(2\)]\) Q5sq \!\(\*SuperscriptBox[\(Q6\), \(2\)]\)-1024 \!\(\*SuperscriptBox[\(Q2\), \(6\)]\) \!\(\*SuperscriptBox[\(Q6\), \(3\)]\)+9216 \!\(\*SuperscriptBox[\(Q2\), \(4\)]\) Q4 \!\(\*SuperscriptBox[\(Q6\), \(3\)]\)-17280 \!\(\*SuperscriptBox[\(Q2\), \(2\)]\) \!\(\*SuperscriptBox[\(Q4\), \(2\)]\) \!\(\*SuperscriptBox[\(Q6\), \(3\)]\)-13824 \!\(\*SuperscriptBox[\(Q4\), \(3\)]\) \!\(\*SuperscriptBox[\(Q6\), \(3\)]\)-32400 Q2 Q5sq \!\(\*SuperscriptBox[\(Q6\), \(3\)]\)-13824 \!\(\*SuperscriptBox[\(Q2\), \(3\)]\) \!\(\*SuperscriptBox[\(Q6\), \(4\)]\)+62208 Q2 Q4 \!\(\*SuperscriptBox[\(Q6\), \(4\)]\)-46656 \!\(\*SuperscriptBox[\(Q6\), \(5\)]\)//Factor"


"Trying to calculate phi(alpha^2)"
phi1 (x1 x2 + x1 x3 + x1 x4 + x1 x5 + x2 x3 + x2 x4 + x2 x5 + x3 x4 + x3 x5 + x4 x5)//Expand;
phi1^2 - (x1^2 x2^2 + x1^2 x3^2 + x1^2 x4^2 + x1^2 x5^2 + x2^2 x3^2 + x2^2 x4^2 + x2^2 x5^2 + x3^2 x4^2 + x3^2 x5^2 + x4^2 x5^2) - 2(x1 x2 x3 x4 + x2 x3 x4 x5 + x3 x4 x5 x1 + x4 x5 x1 x2 + x5 x1 x2 x3)//Expand

"phi(alpha,alpha^2)"
phi1 (x1 + x2 + x3 + x4 + x5) - (x1^2 x2 + x1 x2^2 + x2^2 x3 + x3^2 x2 + x3^2 x4 + x4^2 x3 + x4^2 x5 + x5^2 x4 + x5^2 x1 + x1^2 x5 - x1^2 x3 - x3^2 x1 - x3^2 x5 - x5^2 x3 - x5^2 x2 - x2^2 x5 - x2^2 x4 - x4^2 x2 - x4^2 x1 - x1^2 x4)//Expand;
"SymmetricReduction[(phi1 (x1 + x2 + x3 + x4 + x5) - (x1^2 x2 + x1 x2^2 + x2^2 x3 + x3^2 x2 + x3^2 x4 + x4^2 x3 + x4^2 x5 + x5^2 x4 + x5^2 x1 + x1^2 x5 - x1^2 x3 - x3^2 x1 - x3^2 x5 - x5^2 x3 - x5^2 x2 - x2^2 x5 - x2^2 x4 - x4^2 x2 - x4^2 x1 - x1^2 x4)) sqdiscx,{x1,x2,x3,x4,x5},{s1,s2,s3,s4,s5}]//Expand"



(* ::InheritFromParent:: *)
(**)
