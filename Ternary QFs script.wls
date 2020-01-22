#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]
lambda[i_,j_,k_,l_]:=Det[{{A[[i,j]],B[[i,j]]},{A[[k,l]],B[[k,l]]}}]
A:={{Subscript[a, 11],Subscript[a, 12],Subscript[a, 13]},{Subscript[a, 12],Subscript[a, 22],Subscript[a, 23]},{Subscript[a, 13],Subscript[a, 23],Subscript[a, 33]}}
B:={{Subscript[b, 11],Subscript[b, 12],Subscript[b, 13]},{Subscript[b, 12],Subscript[b, 22],Subscript[b, 23]},{Subscript[b, 13],Subscript[b, 23],Subscript[b, 33]}}

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

m[x_,y_,z_]:=x m1+y m2+z m3//Simplify;
n[x_,y_,z_]:=m[x,y,z]-(Tr[m[x,y,z]]/4 IdentityMatrix[4])//Simplify
tp[x_,y_,z_,p_]:=Tr[MatrixPower[n[x,y,z],p]]//Simplify

"For calculating q in 'TQF ring calc.wls':"
tp[-1,0,0,2]//Expand

m1//MatrixForm
m2//MatrixForm
m3//MatrixForm
"m1 in case of (A0,Bf) conjugated to make alpha the first non-triv basis vector....doesn't work with basis that I thought....:
Note: Now I understand how g in GL3 acts on the basis and see why we get all these minus signs here.  See notes for details."

({
 {0, -q, -s, -r},
 {1, 0, 0, 0},
 {0, 0, 0, 1},
 {0, -1, 0, 0}
})

\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", "q"}], 
RowBox[{"-", "s"}], 
RowBox[{"-", "r"}]},
{"1", "0", "0", "0"},
{"0", "0", "0", 
RowBox[{"-", "1"}]},
{"0", 
RowBox[{"-", "1"}], "0", "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", "s"}], "0", "0"},
{"0", "0", "0", "s"},
{"1", "0", 
RowBox[{"-", "r"}], 
RowBox[{"-", "q"}]},
{"0", "0", "s", "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)
\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", "r"}], "0", 
RowBox[{"-", "s"}]},
{"0", "0", "s", "r"},
{"0", 
RowBox[{"-", "1"}], 
RowBox[{"-", "q"}], "0"},
{"1", "0", "0", 
RowBox[{"-", "q"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)

"t3 = tp[x,y,z,3]//Simplify;
t=Coefficient[t3,x^3]//Simplify"

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
j210[i_,j_,k_]:=Signature[{i,k,j}](Det[{{A1m[[1,i]],A1m[[1,j]],A1m[[2,j]]},{B1m[[1,i]],B1m[[1,j]],B1m[[2,j]]},{cf[[1,i]],cf[[1,j]],cf[[2,j]]}}]+Det[{{A1m[[1,i]],A1m[[2,k]],A1m[[2,i]]},{B1m[[1,i]],B1m[[2,k]],B1m[[2,i]]},{cf[[1,i]],cf[[2,k]],cf[[2,i]]}}])
j111=-Det[{{A1m[[1,1]],A1m[[1,2]],A1m[[1,3]]},{B1m[[1,1]],B1m[[1,2]],B1m[[1,3]]},{cf[[1,1]],cf[[1,2]],cf[[1,3]]}}]-2 Det[{{A1m[[2,1]],A1m[[2,2]],A1m[[2,3]]},{B1m[[2,1]],B1m[[2,2]],B1m[[2,3]]},{cf[[2,1]],cf[[2,2]],cf[[2,3]]}}];
j3[1]

"Calculation of j (corresponds to r; more specifically, r = -8j)"
j[x_,y_,z_]:=j3[1]x^3+j3[2]y^3+j3[3]z^3+j210[1,2,3]x^2 y+j210[1,3,2]x^2 z+j210[2,1,3]y^2 x+j210[2,3,1]y^2 z+j210[3,1,2]z^2 x+j210[3,2,1]z^2 y+j111 x y z
"j[x,y,z]//Expand"
"24 j[x,y,z] - tp[x,y,z,3]//Expand"
"8 j[x,y,z] + Coefficient[CharacteristicPolynomial[n[x,y,z],w],w]//Expand"
Coefficient[CharacteristicPolynomial[n[1,0,0],w],w]

f[x_,y_,z_]:=cf[[1,1]]x^2 + cf[[1,2]]y^2 + cf[[1,3]]z^2 + 2 cf[[2,1]]y z + 2 cf[[2,2]]x z + 2 cf[[2,3]]x y
f[x,y,z]//Expand

aqf[x_,y_,z_]:=A[[1,1]]x^2 + A[[2,2]]y^2 + A[[3,3]]z^2 + A[[1,2]]x y + A[[1,3]]x z + A[[2,3]]y z
bqf[x_,y_,z_]:=B[[1,1]]x^2 + B[[2,2]]y^2 + B[[3,3]]z^2 + B[[1,2]]x y + B[[1,3]]x z + B[[2,3]]y z
cqf[x_,y_,z_]:=(Coefficient[CharacteristicPolynomial[n[x,y,z],w],w^2]+4 (t2 aqf[x,y,z] + t1 bqf[x,y,z]))/3//Expand//Simplify;
fcoeff={{Coefficient[cqf[x,y,z],x^2],Coefficient[cqf[x,y,z],y^2],Coefficient[cqf[x,y,z],z^2]},{0.5 Coefficient[cqf[x,y,z],y z],0.5 Coefficient[cqf[x,y,z],x z],0.5 Coefficient[cqf[x,y,z],x y]}}//Expand;
0.5 cqf[x,y,z] + f[x,y,z] - t2 aqf[x,y,z] - t1 bqf[x,y,z]//Expand;

j[1,0,0]^2 - f[1,0,0]^3 + f[1,0,0]^2 (t1 bqf[1,0,0] + t2 aqf[1,0,0]) - f[1,0,0] (d2 t1 aqf[1,0,0]^2 + d1 t2 bqf[1,0,0]^2 + (t1 t2 - 3 d1 d2) aqf[1,0,0] bqf[1,0,0]) + d2^2 d1 aqf[1,0,0]^3 + d1^2 d2 bqf[1,0,0]^3 - d2 (2 d1 t2 - t1^2) aqf[1,0,0]^2 bqf[1,0,0] - d1 (2 d2 t1 - t2^2) aqf[1,0,0] bqf[1,0,0]^2//Expand

j[x,y,z]^2 - f[x,y,z]^3 + f[x,y,z]^2 (t1 bqf[x,y,z] + t2 aqf[x,y,z]) - f[x,y,z] (d2 t1 aqf[x,y,z]^2 + d1 t2 bqf[x,y,z]^2 + (t1 t2 - 3 d1 d2) aqf[x,y,z] bqf[x,y,z]) + d2^2 d1 aqf[x,y,z]^3 + d1^2 d2 bqf[x,y,z]^3 - d2 (2 d1 t2 - t1^2) aqf[x,y,z]^2 bqf[x,y,z] - d1 (2 d2 t1 - t2^2) aqf[x,y,z] bqf[x,y,z]^2//Expand


"Checking properties hold i.e. code correct"

nphi[x_,y_,z_]:={{fco[x,y,z],-16 (d1 t2 bqf[x,y,z] - d1 d2 aqf[x,y,z]), -16 (t1 d2 aqf[x,y,z] + d1 d2 bqf[x,y,z])},{bqf[x,y,z],fco[x,y,z] - 4 t1 bqf[x,y,z],-4 d2 aqf[x,y,z]},{aqf[x,y,z],-4 d1 bqf[x,y,z],fco[x,y,z]-4 t2 aqf[x,y,z]}}
nphicharcoeff[x_,y_,z_]:={4 t2 aqf[x,y,z] + 4 t1 bqf[x,y,z] - 3 fco[x,y,z],3 fco[x,y,z]^2 - 48 aqf[x,y,z] bqf[x,y,z] d1 d2 - 8 bqf[x,y,z] fco[x,y,z] t1 + 16 aqf[x,y,z]^2 d2 t1 - 8 aqf[x,y,z] fco[x,y,z] t2 + 16 bqf[x,y,z]^2 d1 t2 + 16 aqf[x,y,z] bqf[x,y,z] t1 t2,- fco[x,y,z]^3 + 48 aqf[x,y,z] bqf[x,y,z] fco[x,y,z] d1 d2 + 64 bqf[x,y,z]^3 d1^2 d2 + 64 d1 d2^2 aqf[x,y,z]^3 + 4 bqf[x,y,z] fco[x,y,z]^2 t1 - 16 aqf[x,y,z]^2 fco[x,y,z] d2 t1 - 128 aqf[x,y,z] bqf[x,y,z]^2 d1 d2 t1 + 64 aqf[x,y,z]^2 bqf[x,y,z] d2 t1^2 + 4 aqf[x,y,z] fco[x,y,z]^2 t2 - 16 bqf[x,y,z]^2 fco[x,y,z] d1 t2 - 128 aqf[x,y,z]^2 bqf[x,y,z] d1 d2 t2 - 16 aqf[x,y,z] bqf[x,y,z] fco[x,y,z] t1 t2 + 64 aqf[x,y,z] bqf[x,y,z]^2 d1 t2^2};

nchar[x_,y_,z_]:=CharacteristicPolynomial[n[x,y,z],w]
ncharcoeff[x_,y_,z_]:={Coefficient[nchar[x,y,z],w^2],Coefficient[nchar[x,y,z],w],Coefficient[w nchar[x,y,z],w]}
resolventcoeff[x_,y_,z_]:={-ncharcoeff[x,y,z][[1]],-4 ncharcoeff[x,y,z][[3]],4 ncharcoeff[x,y,z][[1]] ncharcoeff[x,y,z][[3]]-ncharcoeff[x,y,z][[2]]^2}
resolventcoeff[x,y,z]-nphicharcoeff[x,y,z]//Expand;

index1[x_,y_,z_,r_,s_,t_]:=Det[{{1,0,0,0},{0,x,y,z},{0,r,s,t},(m[x,y,z].m[r,s,t])[[All,1]]}]
index2[x_,y_,z_,r_,s_,t_]:=bqf[x,y,z] aqf[r,s,t] - bqf[r,s,t] aqf[x,y,z]
index1[x,y,z,r,s,t] - index2[x,y,z,r,s,t]//Expand

alphap[x_,y_,z_,p_]:=MatrixPower[m[x,y,z],p][[All,1]]
indexq[x_,y_,z_]:=Det[{{1,0,0,0},{0,x,y,z},alphap[x,y,z,2],alphap[x,y,z,3]}]
indexr[x_,y_,z_]:= -4 Det[A1 bqf[x,y,z] - B1 aqf[x,y,z]]
indexq[x,y,z]-indexr[x,y,z]//Expand

discq = Det[{{4,Tr[m1],Tr[m2],Tr[m3]},{Tr[m1],Tr[m1.m1],Tr[m1.m2],Tr[m1.m3]},{Tr[m2],Tr[m1.m2],Tr[m2.m2],Tr[m2.m3]},{Tr[m3],Tr[m3.m1],Tr[m3.m2],Tr[m3.m3]}}];
discr = Discriminant[4 Det[A1 x - B1],x];
discr - discq//Expand


"Computing q and r at [1,0,0] with a_11 = a_12 = a_23 = 0:"
Coefficient[nchar[1,0,0],w^2]//Simplify;
1/8 (-3 \!\(\*SubsuperscriptBox[\(a\), \(13\), \(2\)]\) \!\(\*SubsuperscriptBox[\(b\), \(12\), \(2\)]\)+8 \!\(\*SubsuperscriptBox[\(a\), \(13\), \(2\)]\) Subscript[b, 11] Subscript[b, 22]+8 Subscript[a, 22] (2 Subscript[a, 33] \!\(\*SubsuperscriptBox[\(b\), \(11\), \(2\)]\)-Subscript[a, 13] Subscript[b, 11] Subscript[b, 13]))
Coefficient[nchar[1,0,0],w]//Simplify;
1/8 (-\!\(\*SubsuperscriptBox[\(a\), \(13\), \(3\)]\) (\!\(
\*SubsuperscriptBox[\(b\), \(12\), \(3\)] - \(4\ 
\*SubscriptBox[\(b\), \(11\)]\ 
\*SubscriptBox[\(b\), \(12\)]\ 
\*SubscriptBox[\(b\), \(22\)]\)\))+\!\(
\*SubsuperscriptBox[\(a\), \(13\), \(2\)]\ \((4\ 
\*SubscriptBox[\(a\), \(22\)]\ 
\*SubscriptBox[\(b\), \(11\)]\ \((\(-
\*SubscriptBox[\(b\), \(12\)]\)\ 
\*SubscriptBox[\(b\), \(13\)] + 2\ 
\*SubscriptBox[\(b\), \(11\)]\ 
\*SubscriptBox[\(b\), \(23\)])\))\)\))




(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
