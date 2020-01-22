#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]
lambda[i_,j_,k_,l_]:=Det[{{A[[i,j]],B[[i,j]]},{A[[k,l]],B[[k,l]]}}]
A:={{1,0,0},{0,-a2,a1},{0,a1,a0-a2^2}}
B:={{0,0,1},{0,-1,0},{1,0,-a2}}

A' = X A + r B
B' = Y A + s B 

A1=1/2{{2 Subscript[a, 11],Subscript[a, 12],Subscript[a, 13]},{Subscript[a, 12],2 Subscript[a, 22],Subscript[a, 23]},{Subscript[a, 13],Subscript[a, 23],2 Subscript[a, 33]}};
B1=1/2{{2 Subscript[b, 11],Subscript[b, 12],Subscript[b, 13]},{Subscript[b, 12],2 Subscript[b, 22],Subscript[b, 23]},{Subscript[b, 13],Subscript[b, 23],2 Subscript[b, 33]}};
d1 = Coefficient[Det[A1 x + B1 y],x^3];
d2 = Coefficient[Det[A1 x + B1 y],y^3];
t1 = Coefficient[Det[A1 x + B1 y],x^2 y];
t2 = Coefficient[Det[A1 x + B1 y],x y^2];
pd1={D[d1,Subscript[a, 11]],D[d1,Subscript[a, 22]],D[d1,Subscript[a, 33]],D[d1,Subscript[a, 23]],D[d1,Subscript[a, 13]],D[d1,Subscript[a, 12]]};
pd2={D[d2,Subscript[b, 11]],D[d2,Subscript[b, 22]],D[d2,Subscript[b, 33]],D[d2,Subscript[b, 23]],D[d2,Subscript[b, 13]],D[d2,Subscript[b, 12]]};
cf={{pd1[[2]]pd2[[3]]+pd1[[3]]pd2[[2]]-2 pd1[[4]]pd2[[4]],pd1[[3]]pd2[[1]]+pd1[[1]]pd2[[3]]-2 pd1[[5]]pd2[[5]],pd1[[1]]pd2[[2]]+pd1[[2]]pd2[[1]]-2 pd1[[6]]pd2[[6]]},{pd1[[5]]pd2[[6]]+pd1[[6]]pd2[[5]]-pd1[[1]]pd2[[4]]-pd1[[4]]pd2[[1]],pd1[[6]]pd2[[4]]+pd1[[4]]pd2[[6]]-pd1[[2]]pd2[[5]]-pd1[[5]]pd2[[2]],pd1[[4]]pd2[[5]]+pd1[[5]]pd2[[4]]-pd1[[3]]pd2[[6]]-pd1[[6]]pd2[[3]]}};
A1m={{A1[[1,1]],A1[[2,2]],A1[[3,3]]},{A1[[2,3]],A1[[1,3]],A1[[1,2]]}};
B1m={{B1[[1,1]],B1[[2,2]],B1[[3,3]]},{B1[[2,3]],B1[[1,3]],B1[[1,2]]}};

j3[i_]:=Det[{{A1m[[1,i]],A1m[[2,Mod[i,3]+1]],A1m[[2,Mod[i+1,3]+1]]},{B1m[[1,i]],B1m[[2,Mod[i,3]+1]],B1m[[2,Mod[i+1,3]+1]]},{cf[[1,i]],cf[[2,Mod[i,3]+1]],cf[[2,Mod[i+1,3]+1]]}}]
{Subscript[a,11] -> A'[[1,1]],Subscript[a,12] -> A'[[1,2]],Subscript[a,13] -> A'[[1,3]],Subscript[a,22] -> A'[[2,2]],Subscript[a,23] -> A'[[2,3]],Subscript[a,33] -> A'[[3,3]],Subscript[b,11] -> B'[[1,1]],Subscript[b,12] -> B'[[1,2]],Subscript[b,13] -> B'[[1,3]],Subscript[b,22] -> B'[[2,2]],Subscript[b,23] -> B'[[2,3]],Subscript[b,33] -> B'[[3,3]]}
j3[1]/.%//Expand

ci[i_]:=((-1)^(i+1))lambda[3-Mod[2-i,3],3-Mod[1-i,3],i,i]
ciii[i_,j_,k_]:=Signature[{i,j,k}]lambda[i,k,i,j]+ci[i]
ciij[i_,j_,k_]:=Signature[{i,j,k}]lambda[i,i,i,k]
ciji[i_,j_,k_]:=(Signature[{i,j,k}]lambda[i,k,j,j]+ci[j])/2//Simplify
cijk[i_,j_,k_]:=Signature[{i,j,k}]lambda[j,j,i,i]

cii0[i_,j_,k_]:=ciji[i,k,j]ciij[i,k,j]-ciii[i,k,j]ciji[k,i,j]+cijk[i,k,j]cijk[i,j,k]-ciij[i,j,k]ciji[k,j,i]+ciji[k,i,j]ciji[k,i,j]-ciij[i,k,j]ciii[k,i,j]//Simplify
cij0[i_,j_,k_]:=ciij[j,i,k]ciij[i,j,k]-ciji[i,j,k]ciji[j,i,k]+ciij[j,k,i]cijk[k,i,j]-cijk[i,j,k]ciji[j,k,i]//Simplify

m1=Transpose[{{0,1,0,0},{cii0[1,2,3],ciii[1,2,3],ciij[1,2,3],ciij[1,3,2]},{cij0[1,2,3],ciji[1,2,3],ciji[2,1,3],cijk[1,2,3]},{cij0[1,3,2],ciji[1,3,2],cijk[1,3,2],ciji[3,1,2]}}];
m2=Transpose[{{0,0,1,0},{cij0[1,2,3],ciji[1,2,3],ciji[2,1,3],cijk[1,2,3]},{cii0[2,1,3],ciij[2,1,3],ciii[2,1,3],ciij[2,3,1]},{cij0[2,3,1],cijk[2,3,1],ciji[2,3,1],ciji[3,2,1]}}];
m3=Transpose[{{0,0,0,1},{cij0[1,3,2],ciji[1,3,2],cijk[1,3,2],ciji[3,1,2]},{cij0[2,3,1],cijk[2,3,1],ciji[2,3,1],ciji[3,2,1]},{cii0[3,2,1],ciij[3,1,2],ciij[3,2,1],ciii[3,2,1]}}];

a2 = 2 t2 x + 2 t1 y - 6z;
a2modz = 2 t2 x + 2 t1 y;
a1 = -8w;
a1sq = 64(z^3 - z^2(t1 x + t2 y) + z(d2 t1 x^2 + (t1 t2 - 3 d1 d2)x y + d1 t2 y^2) + (-d2^2 d1 x^3 + d2 (2 d1 t2 - t1^2)x^2 y + d1(2 d2 t1 - t2^2)x y^2 - d1^2 d2 y^3));
a1sqmodz = 64(-d2^2 d1 x^3 + d2 (2 d1 t2 - t1^2)x^2 y + d1(2 d2 t1 - t2^2)x y^2 - d1^2 d2 y^3);
a0 = t2^2 x^2 + t1^2 y^2 + 2 t1 t2 x y + 2 t2 x z + 2 t1 y z - 3z^2 -4(d2 t1 x^2 + (t1 t2-3d1 d2)x y + d1 t2 y^2);
a0modz = t2^2 x^2 + t1^2 y^2 + 2 t1 t2 x y -4(d2 t1 x^2 + (t1 t2-3d1 d2)x y + d1 t2 y^2);
"r = 4 t1 x y - 4 d1 y^2;"

discmodz = -r x(x a0modz - r a2modz)-1/4 x^3 a1sqmodz + 1/4 r (r-x a2modz)^2//Expand;
discmodz - 16 d1 (d1 y^3 - t1 y^2 x + t2 y x^2 - d2 x^3)^2//Expand;











