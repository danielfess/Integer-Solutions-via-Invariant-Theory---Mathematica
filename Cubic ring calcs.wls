#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll["Global`*"]

malpha = Transpose[{{0,1,0},{-16 d1 t2,-4 t1,-4 d1},{16 d1 d2,0,0}}];
mbeta = Transpose[{{0,0,1},{16 d1 d2,0,0},{-16 t1 d2, -4 d2, -4 t2}}];
Collect[-CharacteristicPolynomial[c IdentityMatrix[3] + B malpha + A mbeta,w]//Expand,w]
