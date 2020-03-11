#!/usr/bin/env wolframscript
(* ::Package:: *)

g = {{a^3, a^2 b, a b^2, b^3},{3a^2 c, 2a b c + a^2 d, 2a b d + b^2 c, 3b^2 d},{3a c^2, 2a c d + b c^2, 2b c d + a d^2, 3b d^2},{c^3, c^2 d, c d^2, d^3}};
g//MatrixForm

Q1 = {{0,0,1/2,0},{0,-1,0,0},{1/2,0,0,0},{0,0,0,0}};
Q2 = {{0,0,0,0},{0,0,0,1/2},{0,0,-1,0},{0,1/2,0,0}};
Q3 = {{0,0,0,1/2},{0,0,-1/2,0},{0,-1/2,0,0},{1/2,0,0,0}};
Q4f = {{0,f0/2,0,0},{f0/2,f1,f2/2,0},{0,f2/2,f3,f4/2},{0,0,f4/2,f5}};
Q5f = {{f0,f1/2,0,0},{f1/2,f2,f3/2,0},{0,f3/2,f4,f5/2},{0,0,f5/2,0}};

x1 = a x + c y;
y1 = b x + d y;

G = f0 x1^5 + f1 x1^4 y1 + f2 x1^3 y1^2 + f3 x1^2 y1^3 + f4 x1 y1^4 + f5 y1^5;
g0 = Coefficient[G,x^5]
g1 = Coefficient[G,x^4 y]//Expand
g2 = Coefficient[G,x^3 y^2]//Expand
g3 = Coefficient[G,x^2 y^3]//Expand
g4 = Coefficient[G,x y^4]//Expand
g5 = Coefficient[G,y^5]

Q4g = {{0,g0/2,0,0},{g0/2,g1,g2/2,0},{0,g2/2,g3,g4/2},{0,0,g4/2,g5}};
Q5g = {{g0,g1/2,0,0},{g1/2,g2,g3/2,0},{0,g3/2,g4,g5/2},{0,0,g5/2,0}};

g.Q1.Transpose[g] - (a d - b c)^2 (a^2 Q1 + c^2 Q2 + a c Q3)//Simplify;
g.Q2.Transpose[g] - (a d - b c)^2 (b^2 Q1 + d^2 Q2 + b d Q3)//Simplify;
g.Q3.Transpose[g] - (a d - b c)^2 (2a b Q1 + 2c d Q2 + (b c + a d) Q3)//Simplify;
g.Q4f.Transpose[g] - b Q5g - d Q4g - (3 a^3 b c^2 f0+a c (b c+2 a d) (a^2 f0+2 a b f1+b^2 f2)+b d (2 b c+a d) (a^2 f2+2 a b f3+b^2 f4)+3 b^3 d^2 (a f4+2 b f5))Q1 - (a^3 c^2 d f0+a^2 b c (c^2 f0+2 c d f1+d^2 f2)+a b^2 d (c^2 f2+2 c d f3+d^2 f4)+b^3 d^2 (c f4+2 d f5))Q3 - (a^2 c d (4 c^2 f0+2 c d f1+d^2 f2)+2 a b (c^4 f0+2 c^3 d f1+2 c^2 d^2 f2+2 c d^3 f3+d^4 f4)+b^2 d (c^3 f2+2 c^2 d f3+4 c d^2 f4+6 d^3 f5))Q2//Simplify;

c41 = (3 a^3 b c^2 f0+a c (b c+2 a d) (a^2 f0+2 a b f1+b^2 f2)+b d (2 b c+a d) (a^2 f2+2 a b f3+b^2 f4)+3 b^3 d^2 (a f4+2 b f5))//Simplify;
c43 = (a^3 c^2 d f0+a^2 b c (c^2 f0+2 c d f1+d^2 f2)+a b^2 d (c^2 f2+2 c d f3+d^2 f4)+b^3 d^2 (c f4+2 d f5))//Simplify;
c42 = (a^2 c d (4 c^2 f0+2 c d f1+d^2 f2)+2 a b (c^4 f0+2 c^3 d f1+2 c^2 d^2 f2+2 c d^3 f3+d^4 f4)+b^2 d (c^3 f2+2 c^2 d f3+4 c d^2 f4+6 d^3 f5))//Simplify;

g.Q5f.Transpose[g] - a Q5g - c Q4g - (3 a^3 c^2 (2 a f0+b f1)+a c (b c+2 a d) (a^2 f1+2 a b f2+b^2 f3)+3 a b^3 d^2 f5+b d (2 b c+a d) (a^2 f3+2 a b f4+b^2 f5))Q1 - (a^2 c (6 c^3 f0+4 c^2 d f1+2 c d^2 f2+d^3 f3)+b^2 c d (c^2 f3+2 c d f4+4 d^2 f5)+2 a b (c^4 f1+2 c^3 d f2+2 c^2 d^2 f3+2 c d^3 f4+d^4 f5))Q2 - (a^3 c^2 (2 c f0+d f1)+a^2 b c (c^2 f1+2 c d f2+d^2 f3)+b^3 c d^2 f5+a b^2 d (c^2 f3+2 c d f4+d^2 f5))Q3//Simplify;
c51 = (3 a^3 c^2 (2 a f0+b f1)+a c (b c+2 a d) (a^2 f1+2 a b f2+b^2 f3)+3 a b^3 d^2 f5+b d (2 b c+a d) (a^2 f3+2 a b f4+b^2 f5))//Simplify;
c52 = (a^2 c (6 c^3 f0+4 c^2 d f1+2 c d^2 f2+d^3 f3)+b^2 c d (c^2 f3+2 c d f4+4 d^2 f5)+2 a b (c^4 f1+2 c^3 d f2+2 c^2 d^2 f3+2 c d^3 f4+d^4 f5))//Simplify;
c53 = (a^3 c^2 (2 c f0+d f1)+a^2 b c (c^2 f1+2 c d f2+d^2 f3)+b^3 c d^2 f5+a b^2 d (c^2 f3+2 c d f4+d^2 f5))//Simplify;

A = {{0, -f2 x - f3 y, -f0 w - f1 x, y, z},{f2 x + f3 y, 0, f4 y + f5 z, w, x},{f0 w + f1 x, - f4 y - f5 z, 0, -x, -y},{- y, - w, x, 0, 0},{- z, - x, y, 0, 0}};
subpfaff = Reverse[Diagonal[Minors[A]]]//Factor;

det2 = (a d - b c)^2;
beta = {{a^2 det2, b^2 det2, 2 a b det2, -c41, c51},{c^2 det2, d^2 det2, 2 c d det2, -c42, c52},{a c det2, b d det2, (b c + a d) det2, -c43, c53},{0,0,0,d,-c},{0,0,0,-b,a}};
beta//MatrixForm
columns = beta[[1;;3,4;;5]]//Expand;
columns//MatrixForm

A1f = {{0,0,-f0,0,0},{0,0,0,1,0},{f0,0,0,0,0},{0,-1,0,0,0},{0,0,0,0,0}};
A2f = {{0,-f2,-f1,0,0},{f2,0,0,0,1},{f1,0,0,-1,0},{0,0,1,0,0},{0,-1,0,0,0}};
A3f = {{0,-f3,0,1,0},{f3,0,f4,0,0},{0,-f4,0,0,-1},{-1,0,0,0,0},{0,0,1,0,0}};
A4f = {{0,0,0,0,1},{0,0,f5,0,0},{0,-f5,0,0,0},{0,0,0,0,0},{-1,0,0,0,0}};

Af = {A1f,A2f,A3f,A4f};

A1g = {{0,0,-g0,0,0},{0,0,0,1,0},{g0,0,0,0,0},{0,-1,0,0,0},{0,0,0,0,0}};
A2g = {{0,-g2,-g1,0,0},{g2,0,0,0,1},{g1,0,0,-1,0},{0,0,1,0,0},{0,-1,0,0,0}};
A3g = {{0,-g3,0,1,0},{g3,0,g4,0,0},{0,-g4,0,0,-1},{-1,0,0,0,0},{0,0,1,0,0}};
A4g = {{0,0,0,0,1},{0,0,g5,0,0},{0,-g5,0,0,0},{0,0,0,0,0},{-1,0,0,0,0}};

Ag = {A1g,A2g,A3g,A4g};

Afbeta = {beta.A1f.Transpose[beta],beta.A2f.Transpose[beta],beta.A3f.Transpose[beta],beta.A4f.Transpose[beta]};
Ag - (a d - b c) g.Afbeta//Factor;

"Calculating if there is a form of \Phi equivariant for GL_2:"
columns2 = columns - (a d - b c)^2 {{a^2, b^2, 2 a b},{c^2, d^2, 2 c d},{a c, b d, b c + a d}}.{{2 f1/5, -4 f2/10},{4 f3/10, -2 f4/5},{f2/10, -f3/10}} + {{2 g1/5, -4 g2/10},{4 g3/10, -2 g4/5},{g2/10, -g3/10}}.{{d,-c},{-b,a}}//Expand;
columns2//MatrixForm

"betaf alters \Phi(f) = Af into Afnew"
betaf = {{1,0,0,2 f1/5, -4 f2/10},{0,1,0,4 f3/10, -2 f4/5},{0,0,1,f2/10, -f3/10},{0,0,0,1,0},{0,0,0,0,1}};
betaf//MatrixForm
A1fnew = betaf.A1f.Transpose[betaf];
A2fnew = betaf.A2f.Transpose[betaf];
A3fnew = betaf.A3f.Transpose[betaf];
A4fnew = betaf.A4f.Transpose[betaf];
A1fnew//MatrixForm
A2fnew//MatrixForm
A3fnew//MatrixForm
A4fnew//MatrixForm
Afnew = {A1fnew,A2fnew,A3fnew,A4fnew};

"betag alters \Phi(g) = Ag into Agnew"
betag = {{1,0,0,2 g1/5, -4 g2/10},{0,1,0,4 g3/10, -2 g4/5},{0,0,1,g2/10, -g3/10},{0,0,0,1,0},{0,0,0,0,1}};
A1gnew = betag.A1g.Transpose[betag];
A2gnew = betag.A2g.Transpose[betag];
A3gnew = betag.A3g.Transpose[betag];
A4gnew = betag.A4g.Transpose[betag];
Agnew = {A1gnew,A2gnew,A3gnew,A4gnew};
Agnew[[1]][[2,4]]

"betahom is the map from GL2 to GL5 for which we hope the new \Phi map is equivariant"
betahom = {{a^2 det2, b^2 det2, 2 a b det2, 0,0},{c^2 det2, d^2 det2, 2 c d det2, 0,0},{a c det2, b d det2, (b c + a d) det2, 0,0},{0,0,0,d,-c},{0,0,0,-b,a}};

"Image of Afnew under betahom"
A1fnewbeta = betahom.A1fnew.Transpose[betahom];
A2fnewbeta = betahom.A2fnew.Transpose[betahom];
A3fnewbeta = betahom.A3fnew.Transpose[betahom];
A4fnewbeta = betahom.A4fnew.Transpose[betahom];

Afbetahom = {A1fnewbeta,A2fnewbeta,A3fnewbeta,A4fnewbeta};
Afbetahom2 = {betahom.A1fnew.Transpose[betahom], betahom.A2fnew.Transpose[betahom], betahom.A3fnew.Transpose[betahom], betahom.A4fnew.Transpose[betahom]};
Afbetahom2 - Afbetahom//Simplify

"Testing if Agnew = image of Afnew under GL2 action via map to GL4 x GL5"
Agnew - (a d - b c) g.Afbetahom//Factor







