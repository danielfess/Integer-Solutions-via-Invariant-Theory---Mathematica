#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]

m = {{0,0,0,-f4/f0},{1,0,0,-f3/f0},{0,1,0,-f2/f0},{0,0,1,-f1/f0}};
m2 = MatrixPower[m,2];
m3 = MatrixPower[m,3];
m4 = MatrixPower[m,4];
m5 = MatrixPower[m,5];
m6 = MatrixPower[m,6];
M = {{4,Tr[m],Tr[m2],Tr[m3]},{Tr[m],Tr[m2],Tr[m3],Tr[m4]},{Tr[m2],Tr[m3],Tr[m4],Tr[m5]},{Tr[m3],Tr[m4],Tr[m5],Tr[m6]}};
P = {{1,0,0,0},{0,f0,f1,f2},{0,0,f0,f1},{0,0,0,f0}};
T = Transpose[P].M.P //Simplify
CofactorT = Inverse[T] Det[T]//Expand;
Dual = CofactorT.Transpose[P]//Expand;
m.Transpose[Dual][[All,4]] - Transpose[Dual][[All,3]]//Expand
m.Transpose[Dual][[All,3]] - Transpose[Dual][[All,2]]//Expand



