#!/usr/bin/env wolframscript
(* ::Package:: *)

A1 = {{0,-((2 f1)/5),-f0,0,0},{(2 f1)/5,0,f2/10,1,0},{f0,-(f2/10),0,0,0},{0,-1,0,0,0},{0,0,0,0,0}};
A2 = {{0,-((3 f2)/5),-((3 f1)/5),0,0},{(3 f2)/5,0,(3 f3)/10,0,1},{(3 f1)/5,-((3 f3)/10),0,-1,0},{0,0,1,0,0},{0,-1,0,0,0}};
A3 = {{0,-((3 f3)/5),-((3 f2)/10),1,0},{(3 f3)/5,0,(3 f4)/5,0,0},{(3 f2)/10,-((3 f4)/5),0,0,-1},{-1,0,0,0,0},{0,0,1,0,0}};
A4 = {{0,-((2 f4)/5),-(f3/10),0,1},{(2 f4)/5,0,f5,0,0},{f3/10,-f5,0,0,0},{0,0,0,0,0},{-1,0,0,0,0}};

invt[k1_,l1_,k2_,l2_,k3_,l3_,k4_,l4_] := Det[{{A1[[k1,l1]], A2[[k1,l1]], A3[[k1,l1]], A4[[k1,l1]]}, {A1[[k2,l2]],A2[[k2,l2]],A3[[k2,l2]],A4[[k2,l2]]},{A1[[k3,l3]],A2[[k3,l3]],A3[[k3,l3]],A4[[k3,l3]]},{A1[[k4,l4]],A2[[k4,l4]],A3[[k4,l4]],A4[[k4,l4]]}}]

invt[1,5,1,4,2,1,2,3];

signproduct[i_,j_,k_,l10_,k20_,l20_,k30_,l30_,k40_,l40_,k11_,l11_,k21_,l21_,k31_,l31_,k41_,l41_,k12_,l12_,k22_,l22_,k32_,l32_,k42_,l42_] :=  Signature[{i,l11,l21,k31,l31}] * Signature[{j,l12,l22,k32,l32}] * Signature[{k11,l10,l20,k30,l30}] * Signature[{k12,k21,k41,l41,l42}] * Signature[{k20,k22,k40,k42,l40}]
invtproduct[i_,j_,k_,l10_,k20_,l20_,k30_,l30_,k40_,l40_,k11_,l11_,k21_,l21_,k31_,l31_,k41_,l41_,k12_,l12_,k22_,l22_,k32_,l32_,k42_,l42_] := invt[k,l10,k20,l20,k30,l30,k40,l40] * invt[k11,l11,k21,l21,k31,l31,k41,l41] * invt[k12,l12,k22,l22,k32,l32,k42,l42]
term[i_,j_,k_,l10_,k20_,l20_,k30_,l30_,k40_,l40_,k11_,l11_,k21_,l21_,k31_,l31_,k41_,l41_,k12_,l12_,k22_,l22_,k32_,l32_,k42_,l42_] := signproduct[i,j,k,l10,k20,l20,k30,l30,k40,l40,k11,l11,k21,l21,k31,l31,k41,l41,k12,l12,k22,l22,k32,l32,k42,l42] * invtproduct[i,j,k,l10,k20,l20,k30,l30,k40,l40,k11,l11,k21,l21,k31,l31,k41,l41,k12,l12,k22,l22,k32,l32,k42,l42]

signproduct[1,1,3,2,1,3,4,5,3,5,1,2,2,3,4,5,3,4,1,2,4,3,4,5,2,5]

"signproduct[i,j,k,l10,k20,l20,k30,l30,k40,l40,k11,l11,k21,l21,k31,l31,k41,l41,k12,l12,k22,l22,k32,l32,k42,l42]"

mult[i_,j_,k_] := 1/2304 * Sum[term[i,j,k,l10,k20,l20,k30,l30,k40,l40,k11,l11,k21,l21,k31,l31,k41,l41,k12,l12,k22,l22,k32,l32,k42,l42],{l10,1,5},{k20,1,5},{l20,1,5},{k30,1,5},{l30,1,5},{k40,1,5},{l40,1,5},{k11,1,5},{l11,1,5},{k21,1,5},{l21,1,5},{k31,1,5},{l31,1,5},{k41,1,5},{l41,1,5},{k12,1,5},{l12,1,5},{k22,1,5},{l22,1,5},{k32,1,5},{l32,1,5},{k42,1,5},{l42,1,5}]
"mult is INCREDIBLY slow.  Unusably slow.  Need to figure out some shortcuts."

"Cubic in fi:"
invt[1,2,1,3,2,3,2,4]
invt[1,2,1,3,2,3,2,5]
invt[1,2,1,3,2,3,1,4]
invt[1,2,1,3,2,3,1,5]

"Quadratic in fi (some also with indices [3,4] and [3,5], which come out identical as [2,5] and [1,4] terms respectively):"
invt[1,2,1,3,2,4,2,5]
invt[1,2,1,3,2,4,1,4]
invt[1,2,1,3,2,4,1,5]
invt[1,2,1,3,2,5,1,4]
invt[1,2,1,3,2,5,1,5]
invt[1,2,1,3,1,4,1,5]
invt[1,2,2,3,2,4,2,5]
invt[1,2,2,3,2,4,1,4]
invt[1,2,2,3,2,4,1,5]
invt[1,2,2,3,2,5,1,4]
invt[1,2,2,3,2,5,1,5]
invt[1,2,2,3,1,4,1,5]
invt[2,3,1,3,2,4,2,5]
invt[2,3,1,3,2,4,1,4]
invt[2,3,1,3,2,4,1,5]
invt[2,3,1,3,2,5,1,4]
invt[2,3,1,3,2,5,1,5]
invt[2,3,1,3,1,4,1,5]

"Linear in fi:"
invt[1,2,2,4,2,5,1,4]
invt[1,2,2,4,2,5,1,5]
invt[1,2,2,4,1,4,1,5]
invt[1,2,2,5,1,4,1,5]
invt[1,3,2,4,2,5,1,4]
invt[1,3,2,4,2,5,1,5]
invt[1,3,2,4,1,4,1,5]
invt[1,3,2,5,1,4,1,5]
invt[2,3,2,4,2,5,1,4]
invt[2,3,2,4,2,5,1,5]
invt[2,3,2,4,1,4,1,5]
invt[2,3,2,5,1,4,1,5]

"Constant i.e. degree zero in fi:"
invt[2,4,2,5,1,4,1,5]





