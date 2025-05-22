#!/usr/bin/wolframscript wolframscript
(* ::Package:: *)

(*SetDirectory[NotebookDirectory[]]*)
LamA=Import["Lamst4.wdx"];
Print["load OK"]

LamB=LamA*quB*qdB
Lam=Table[FullSimplify[LamB[[i]]],{i,10}];


FortranForm[{
dtLAM[1]==Lam[[1]],
dtLAM[2]==Lam[[2]],
dtLAM[3]==Lam[[3]],
dtLAM[4]==Lam[[4]],
dtLAM[5]==Lam[[5]],
dtLAM[6]==Lam[[6]],
dtLAM[7]==Lam[[7]],
dtLAM[8]==Lam[[8]],
dtLAM[9]==Lam[[9]],
dtLAM[10]==Lam[[10]]
}]>>"flow_4quarkst5.f90";
Export["Lamst5.wdx",Lam,"WDX"];


rep4pointout={lambdaVmA->lam[1],
lambdaVpA->lam[2],
lambdaVmAadj->lam[3],
lambdaSigma->lam[4],
lambdaPion->lam[5],
lambdaEta->lam[6],
lambdaA->lam[7],
lambdaSmPmadj->lam[8],
lambdaSpPmadj->lam[9],
lambdaSpPpadj->lam[10]};
