#!/usr/bin/env wolframscript
(* ::Package:: *)

"Need to compute structure coeffs and print mult tables";
q11 = {{0,0,0,0,2 f0}};
q22 = {{-2,0,0,-2 f1, 2 f2}};
q33 = {{0,-2,0, -2 f3, 2 f4}};
q44 = {{0,0,0, -2 f5, 0 }};
q12 = {{0,0,0,-f0,f1}};
q13 = {{1,0,0,0,0}};
q14 = {{0,0,1,0,0}};
q23 = {{0,0,-1,-f2,f3}};
q24 = {{0,1,0,0,0}};
q34 = {{0,0,0,-f4,f5}};

q = {{q11,q12,q13,q14},{q12,q22,q23,q24},{q13,q23,q33,q34},{q14,q24,q34,q44}};

A1 = {{0,0,-f0,0,0},{0,0,0,1,0},{f0,0,0,0,0},{0,-1,0,0,0},{0,0,0,0,0}};
A2 = {{0,-f2,-f1,0,0},{f2,0,0,0,1},{f1,0,0,-1,0},{0,0,1,0,0},{0,-1,0,0,0}};
A3 = {{0,-f3,0,1,0},{f3,0,f4,0,0},{0,-f4,0,0,-1},{-1,0,0,0,0},{0,0,1,0,0}};
A4 = {{0,0,0,0,1},{0,0,f5,0,0},{0,-f5,0,0,0},{0,0,0,0,0},{-1,0,0,0,0}};

A = {A1,A2,A3,A4};

invariant[i_,j_,k_,l_,m_] := (q[[i,j]].A[[k]].Transpose[q[[l,m]]])[[1,1]]

cijk[i_,j_,k_,l_] := Signature[{i,j,k,l}] invariant[i,i,l,j,j]/4
ciij[i_,j_,k_,l_] := Signature[{i,j,k,l}] invariant[l,i,i,i,k]

"Testing cijk:"
cijk[2,3,4,1]
cijk[3,2,4,1]

"For cijj to give the right answer, we need to pick k,l correctly given i,j
i.e. with c_ik^k = 0, and l the remaining index.
Hence we devise the following, so that the (i,j)-th entry gives the needed (i,j,k,l).
Note:  When (i,j) = (i,i) or (i,k), the formula for cijj will spit out 0."
ones = {1,1,1,1};
indices = {{ones,ones,{1,3,2,4},{1,4,2,3}},{ones,ones,{2,3,1,4},{2,4,1,3}},{{3,1,4,2},{3,2,4,1},ones,ones},{{4,1,3,2},{4,2,3,1},ones,ones}};
cijj[i_,j_] := Signature[indices[[i,j]]] invariant[j,indices[[i,j]][[3]],indices[[i,j]][[4]],i,i]/2

"Testing cijj and indices.  cijj works even for bad (i,j) since the Signature comes out as 0."
indices[[1,3]]
cijj[1,2]
cijj[3,1]
cijj[3,2]

"To calculate c_ii^i we also introduce fixed indices for each i.
k is fixed given i, so there are two choices of indices for each i:"
indices2 = {{1,3,2,4},{2,3,1,4},{3,1,4,2},{4,1,3,2}}
indices3 = {{1,4,2,3},{2,4,1,3},{3,2,4,1},{4,2,3,1}}

ciii2[i_] := Signature[indices2[[i]]](invariant[indices2[[i]][[2]],indices2[[i]][[3]],indices2[[i]][[4]],i,i]/2 + invariant[i,indices2[[i]][[2]],indices2[[i]][[4]],indices2[[i]][[3]],i])
ciii3[i_] := Signature[indices3[[i]]](invariant[indices3[[i]][[2]],indices3[[i]][[3]],indices3[[i]][[4]],i,i]/2 + invariant[i,indices3[[i]][[2]],indices3[[i]][[4]],indices3[[i]][[3]],i])

"Tests below indicate ciii2 and ciii3 give the same answers."
ciii2[1]
ciii3[1]
ciii2[2]
ciii3[2]
ciii2[3]
ciii3[3]
ciii2[4]
ciii3[4]

"cij0 next. Need to separate into i=j and i!=j cases?
ind4ii has (i,k,..).  ind4ij has (i,j,k,..)."
ind4ii = {{1,2,3,4},{2,3,4,1},{3,4,1,2},{4,1,2,3}}
ind4ij = {{ones,{1,2,3,4},{1,3,2,4},{1,4,2,3}},{{}}}

cii0[i_] := cijj[Mod[i,4]+1,i] ciij[i,Mod[i,4]+1,Mod[i+1,4]+1,Mod[i+2,4]+1] - ciii2[i] cijj[i,Mod[i,4]+1] + cijj[i,Mod[i,4]+1]^2 - ciij[i,Mod[i,4]+1,Mod[i+1,4]+1,Mod[i+2,4]+1] ciii2[Mod[i,4]+1] + cijk[i,Mod[i,4]+1,Mod[i+1,4]+1,Mod[i+2,4]+1] cijk[i,Mod[i+1,4]+1,Mod[i,4]+1,Mod[i+2,4]+1] - ciij[i,Mod[i+1,4]+1,Mod[i,4]+1,Mod[i+2,4]+1] cijj[Mod[i+1,4]+1,Mod[i,4]+1] + cijk[i,Mod[i,4]+1,Mod[i+2,4]+1,Mod[i+1,4]+1] cijk[i,Mod[i+2,4]+1,Mod[i,4]+1,Mod[i+1,4]+1] - ciij[i,Mod[i+2,4]+1,Mod[i+1,4]+1,Mod[i,4]+1] cijj[Mod[i+2,4]+1,Mod[i,4]+1]

cii0[1]
cii0[2]
cii0[3]
cii0[4]

"Indices for r when i !=j :"
ind4ij = {{{0,0},{3,4},{2,4},{2,3}},{{3,4},{0,0},{1,4},{1,3}},{{2,4},{1,4},{0,0},{1,2}},{{2,3},{1,3},{1,2},{0,0}}}

cij0[i_,j_] := ciij[j,i,ind4ij[[i,j]][[1]],ind4ij[[i,j]][[2]]] ciij[i,j,ind4ij[[i,j]][[1]],ind4ij[[i,j]][[2]]] - cijj[j,i] cijj[i,j] + ciij[j,ind4ij[[i,j]][[1]],ind4ij[[i,j]][[2]],i] cijk[i,ind4ij[[i,j]][[1]],j,ind4ij[[i,j]][[2]]] - cijk[i,j,ind4ij[[i,j]][[1]],ind4ij[[i,j]][[2]]] cijj[ind4ij[[i,j]][[1]],j] + ciij[j,ind4ij[[i,j]][[2]],ind4ij[[i,j]][[1]],i] cijk[i,ind4ij[[i,j]][[2]],j,ind4ij[[i,j]][[1]]] - cijk[i,j,ind4ij[[i,j]][[2]],ind4ij[[i,j]][[1]]] cijj[ind4ij[[i,j]][[2]],j]

cij0[1,2]
cij0[2,1]
cij0[1,3]
cij0[3,1]

m1 = Transpose[{{0,1,0,0,0},{cii0[1],ciii2[1],ciij[1,2,3,4],ciij[1,3,2,4],ciij[1,4,2,3]},{cij0[1,2],cijj[2,1],cijj[1,2],cijk[1,2,3,4],cijk[1,2,4,3]},{cij0[1,3],cijj[3,1],cijk[1,3,2,4],cijj[1,3],cijk[1,3,4,2]},{cij0[1,4],cijj[4,1],cijk[1,4,2,3],cijk[1,4,3,2],cijj[1,4]}}];
m2 = Transpose[{{0,0,1,0,0},{cij0[1,2],cijj[2,1],cijj[1,2],cijk[1,2,3,4],cijk[1,2,4,3]},{cii0[2],ciij[2,1,3,4],ciii2[2],ciij[2,3,1,4],ciij[2,4,1,3]},{cij0[2,3],cijk[2,3,1,4],cijj[3,2],cijj[2,3],cijk[2,3,4,1]},{cij0[2,4],cijk[2,4,1,3],cijj[4,2],cijk[2,4,3,1],cijj[2,4]}}];
m3 = Transpose[{{0,0,0,1,0},{cij0[1,3],cijj[3,1],cijk[1,3,2,4],cijj[1,3],cijk[1,3,4,2]},{cij0[2,3],cijk[2,3,1,4],cijj[3,2],cijj[2,3],cijk[2,3,4,1]},{cii0[3],ciij[3,1,2,4],ciij[3,2,1,4],ciii2[3],ciij[3,4,1,2]},{cij0[3,4],cijk[3,4,1,2],cijk[3,4,2,1],cijj[4,3],cijj[3,4]}}];
m4 = Transpose[{{0,0,0,0,1},{cij0[1,4],cijj[4,1],cijk[1,4,2,3],cijk[1,4,3,2],cijj[1,4]},{cij0[2,4],cijk[2,4,1,3],cijj[4,2],cijk[2,4,3,1],cijj[2,4]},{cij0[3,4],cijk[3,4,1,2],cijk[3,4,2,1],cijj[4,3],cijj[3,4]},{cii0[4],ciij[4,1,2,3],ciij[4,2,1,3],ciij[4,3,1,2],ciii2[4]}}];

m1//MatrixForm
m2//MatrixForm
m3//MatrixForm
m4//MatrixForm

"Check matrices multiply commutatively:"
m1.m2 - m2.m1//Simplify;
m1.m3 - m3.m1//Simplify;
m1.m4 - m4.m1//Simplify;
m2.m3 - m3.m2//Simplify;
m2.m4 - m4.m2//Simplify;
m3.m4 - m4.m3//Simplify;

"Check associativity holds:"
m = {m1,m2,m3,m4};
products = {{m1.m1,m1.m2,m1.m3,m1.m4},{m2.m1,m2.m2,m2.m3,m2.m4},{m3.m1,m3.m2,m3.m3,m3.m4},{m4.m1,m4.m2,m4.m3,m4.m4}};
assoc[i_,j_,k_] := m[[i]].products[[j,k]] - products[[i,j]].m[[k]]//Simplify

assoc[1,1,1];
assoc[1,1,2];
assoc[1,1,3];
assoc[1,1,4];
assoc[1,3,1];
assoc[2,4,3];
assoc[4,1,2];
assoc[3,1,3];

"Looking at traces tells us that if this is Q_f then to obtain Nakagawa's basis, we need to translate eta_2 \mapsto eta_2 - f2, and eta_4 \mapsto eta_4 - f4:"
P = {{1,0,f2,0,f4},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}};
m0 = IdentityMatrix[5];
M1 = P.m1.Inverse[P]//MatrixForm
M2 = P.(m2- f2 m0) .Inverse[P]//MatrixForm
M3 = P.m3.Inverse[P]//MatrixForm
M4 = P.(m4 - f4 m0).Inverse[P]//MatrixForm



