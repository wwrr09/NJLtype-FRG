#!/usr/bin/wolframscript wolframscript
(* ::Package:: *)

(* ::Chapter:: *)
(*NJL AMM*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)


(* ::Section:: *)
(*FormTracer*)


PacletManagerAllowInternetInitialValue = PacletManager`$AllowInternet;
PacletManager`$AllowInternet = False;
Needs["FormTracer`"]
PacletManager`$AllowInternet = PacletManagerAllowInternetInitialValue;


DefineLorentzTensors[
delta[lor,a,b](*Lorentz delta*),
vec[p,\[Mu]](*vector*),
sp[p,q](*inner product*),
eps[](*Levi-Civita symbol*),
delta[Spinor,i,j](*Kronecker delta in spinor space*),
gamma[\[Mu],i,j](*Gamma matrices*),
gamma5[i,j](*fifth gamma matrix*),
vecs[p,mu](*spacial vector*),
sps[p,q](*spacial vector product*)
];
DefineGroupTensors[{
{SUNfund,{color(*,3*)},delta[adj,a,b],SUNF[a,b,c],delta[FundCol,a,b],SUNT[a,i,j]},
{SUNfund,{flavor(*,2*)},delta[MesonFlav,a,b],SUNFF[a,b,c],delta[FundFlav,a,b],SUNTF[a,i,j]}
}];


DefineFormAutoDeclareFunctions[{lambda,Slambda,Alambda,p,q,Z,m,tilZ,dtZ,R,dR,qfB,k,tanh,s}];
AddExtraVars[ZQG1,ZQG2,ZQG3,ZQG4,ZQG5,ZQG6,ZQG7,ZQG8];
AddExtraVars[ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj];


(* ::Section:: *)
(*Inner Product*)


repvec={vec[a_,b_]:>ToExpression[ToString[a]<>"t[["<>ToString[b+1]<>"]]"]};
pt={p,0,0,0};
\[Theta]1=(2\[Pi])/3;
\[Theta]2=(4\[Pi])/3;
p1t={p Cos[\[Theta]2],0,0,p Sin[\[Theta]2]};
p2t={p Cos[\[Theta]1],0,0,p Sin[\[Theta]1]};
p3t={p,0,0,0};
qt={q0,q1,q2,q3};
momentumProductionResult={
sp[q,q]->q^2,
sp[a_,b_]:>ToExpression[ToString[a]<>"t."<>ToString[b]<>"t"]
};
innerProductRepl[expr_]:=expr//.{sp[a_+b_,c_]->sp[a,c]+sp[b,c]}//.{sp[a_,b_+c_]->sp[a,b]+sp[a,c]}/.{sp[c_Integer  a_,b_]->c sp[a,b]}/.{sp[a_,c_Integer  b_]->c sp[a,b]}//.{vec[a_+b_,c_]->vec[a,c]+vec[b,c]}/.{vec[c_Integer  a_,b_]->c vec[a,b]}/.repvec/.momentumProductionResult

repvec4={vec[a_,b_]:>ToExpression[ToString[a]<>"t4Ave[["<>ToString[b+1]<>"]]"]};
pt4Ave={p,0,0,0};
p1t4Ave=-pt4Ave;
p2t4Ave=pt4Ave;
p3t4Ave=pt4Ave;
p4t4Ave=-pt4Ave;
qt4Ave={q x1,q Sqrt[1-x1^2]x2,q Sqrt[1-x1^2] Sqrt[1-x2^2]x3,q Sqrt[1-x1^2] Sqrt[1-x2^2] Sqrt[1-x3^2]};
momentumProductionResult4Ave={
sp[q,q]->q^2,
sp[a_,b_]:>ToExpression[ToString[a]<>"t4Ave."<>ToString[b]<>"t4Ave"]
};
innerProductRepl4Ave[expr_]:=expr//.{sp[a_+b_,c_]->sp[a,c]+sp[b,c]}//.{sp[a_,b_+c_]->sp[a,b]+sp[a,c]}/.{sp[c_Integer  a_,b_]->c sp[a,b]}/.{sp[a_,c_Integer  b_]->c sp[a,b]}//.{vec[a_+b_,c_]->vec[a,c]+vec[b,c]}/.{vec[c_Integer  a_,b_]->c vec[a,b]}/.repvec4/.momentumProductionResult4Ave


(* ::Section:: *)
(*Feynman Rules*)


PiPerp[\[Mu]_,\[Nu]_,p_]:=delta[lor,\[Mu],\[Nu]]-(vec[p,\[Mu]]vec[p,\[Nu]])/sp[p,p]


Gq[i_,j_,q_,qf_,s_]:=qexp[qf,s]*(-I (1-tanh[qf,s]^2) (1+Rq[q,s,qf]) (gamma[1,i,j] vec[q,1]+gamma[2,i,j] vec[q,2])+(delta[Spinor,Si,j]+I tanh[qf,s] gamma[1,Si,Sj] gamma[2,Sj,j]) (delta[Spinor,i,Si] mq[qf]-I (1+Rq[q,s,qf]) (gamma[0,i,Si] vec[q,0]+gamma[3,i,Si] vec[q,3])))/.{Si->Unique[Si],Sj->Unique[Sj]};
dGq[i_,j_,q_,qf_,s_]:=qexp[qf,s]*dRq[q,s,qf]*(-I (-1+tanh[qf,s]^2) (gamma[1,i,j] vec[q,1]+gamma[2,i,j] vec[q,2])+I (delta[Spinor,Si,j]+I tanh[qf,s] gamma[1,Si,Sj] gamma[2,Sj,j]) (gamma[0,i,Si] vec[q,0]+gamma[3,i,Si] vec[q,3])+2 s (1+Rq[q,s,qf]) (vec[q,0]^2+(tanh[qf,s] (vec[q,1]^2+vec[q,2]^2))/(qf*s)+vec[q,3]^2) (I (-1+tanh[qf,s]^2) (1+Rq[q,s,qf]) (gamma[1,i,j] vec[q,1]+gamma[2,i,j] vec[q,2])+(delta[Spinor,Si,j]+I tanh[qf,s] gamma[1,Si,Sj] gamma[2,Sj,j]) (delta[Spinor,i,Si] mq[qf]-I (1+Rq[q,s,qf]) (gamma[0,i,Si] vec[q,0]+gamma[3,i,Si] vec[q,3]))))/.{Si->Unique[Si],Sj->Unique[Sj]};


Gqfull[Fcb1_,Fca1_,Fg1_,Ff1_,j1_,i1_,q_,s_]:=delta[FundCol,Fcb1,Fca1]*(delta[FundFlav,Fg1,1]Gq[j1,i1,q,qu,s]delta[FundFlav,1,Ff1]+delta[FundFlav,Fg1,2]Gq[j1,i1,q,qd,s]delta[FundFlav,2,Ff1]);
dGqfull[Fcb1_,Fca1_,Fg1_,Ff1_,j1_,i1_,q_,s_]:=delta[FundCol,Fcb1,Fca1]*(delta[FundFlav,Fg1,1]dGq[j1,i1,q,qu,s]delta[FundFlav,1,Ff1]+delta[FundFlav,Fg1,2]dGq[j1,i1,q,qd,s]delta[FundFlav,2,Ff1]);


Gammaq4[FCa_,Ff_,Si_,FCb_,Fg_,Sj_,FCc_,Fh_,Sk_,FCd_,Fi_,Sl_,p1_,p2_,p3_,p4_]:=-(1/4 ((((((-2 delta[FundCol,FCa,FCd]) delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,Sl]) gamma[SQCDmu,Sk,Sj]+(((((2 delta[FundCol,FCa,FCb]) delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,Sj]) gamma[SQCDmu,Sk,Sl]+((((((delta[FundCol,FCa,FCb] delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,SQCDSn]) gamma[SQCDmu,Sk,SQCDSm]) gamma5[SQCDSm,Sl]) gamma5[SQCDSn,Sj]-((((((delta[FundCol,FCa,FCd] delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,SQCDSm]) gamma[SQCDmu,Sk,SQCDSn]) gamma5[SQCDSm,Sl]) gamma5[SQCDSn,Sj]-((((((delta[FundCol,FCa,FCd] delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,SQCDSn]) gamma[SQCDmu,Sk,SQCDSm]) gamma5[SQCDSm,Sj]) gamma5[SQCDSn,Sl]+((((((delta[FundCol,FCa,FCb] delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,SQCDSm]) gamma[SQCDmu,Sk,SQCDSn]) gamma5[SQCDSm,Sj]) gamma5[SQCDSn,Sl]) lambdaVmA[p1,p2,p3,p4]
+1/4 ((((((-2 delta[FundCol,FCa,FCd]) delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,Sl]) gamma[SQCDmu,Sk,Sj]+(((((2 delta[FundCol,FCa,FCb]) delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,Sj]) gamma[SQCDmu,Sk,Sl]-((((((delta[FundCol,FCa,FCb] delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,SQCDSn]) gamma[SQCDmu,Sk,SQCDSm]) gamma5[SQCDSm,Sl]) gamma5[SQCDSn,Sj]+((((((delta[FundCol,FCa,FCd] delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,SQCDSm]) gamma[SQCDmu,Sk,SQCDSn]) gamma5[SQCDSm,Sl]) gamma5[SQCDSn,Sj]+((((((delta[FundCol,FCa,FCd] delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,SQCDSn]) gamma[SQCDmu,Sk,SQCDSm]) gamma5[SQCDSm,Sj]) gamma5[SQCDSn,Sl]-((((((delta[FundCol,FCa,FCb] delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,SQCDSm]) gamma[SQCDmu,Sk,SQCDSn]) gamma5[SQCDSm,Sj]) gamma5[SQCDSn,Sl]) lambdaVpA[p1,p2,p3,p4]
+1/4 lambdaVmAadj[p1,p2,p3,p4] ((((((-2 delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,Sl]) gamma[SQCDmu,Sk,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]-((((((delta[FundFlav,Ff,Fi] delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,SQCDSm]) gamma[SQCDmu,Sk,SQCDSn]) gamma5[SQCDSm,Sl]) gamma5[SQCDSn,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]-((((((delta[FundFlav,Ff,Fi] delta[FundFlav,Fg,Fh]) gamma[SQCDmu,Si,SQCDSn]) gamma[SQCDmu,Sk,SQCDSm]) gamma5[SQCDSm,Sj]) gamma5[SQCDSn,Sl]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]+(((((2 delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,Sj]) gamma[SQCDmu,Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]+((((((delta[FundFlav,Ff,Fg] delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,SQCDSn]) gamma[SQCDmu,Sk,SQCDSm]) gamma5[SQCDSm,Sl]) gamma5[SQCDSn,Sj]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]+((((((delta[FundFlav,Ff,Fg] delta[FundFlav,Fh,Fi]) gamma[SQCDmu,Si,SQCDSm]) gamma[SQCDmu,Sk,SQCDSn]) gamma5[SQCDSm,Sj]) gamma5[SQCDSn,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd])
+1/4 ((((((-2 delta[FundCol,FCa,FCd]) delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]+(((((2 delta[FundCol,FCa,FCb]) delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) lambdaSigma[p1,p2,p3,p4]
+1/4 lambdaPion[p1,p2,p3,p4] ((((((8 delta[FundCol,FCa,FCd]) delta[FundCol,FCb,FCc]) gamma5[Si,Sl]) gamma5[Sk,Sj]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]-(((((8 delta[FundCol,FCa,FCb]) delta[FundCol,FCc,FCd]) gamma5[Si,Sj]) gamma5[Sk,Sl]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi])
+1/4 ((((((2 delta[FundCol,FCa,FCd]) delta[FundCol,FCb,FCc]) delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma5[Si,Sl]) gamma5[Sk,Sj]-(((((2 delta[FundCol,FCa,FCb]) delta[FundCol,FCc,FCd]) delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma5[Si,Sj]) gamma5[Sk,Sl]) lambdaEta[p1,p2,p3,p4]
+1/4 lambdaA[p1,p2,p3,p4] ((((((-8 delta[FundCol,FCa,FCd]) delta[FundCol,FCb,FCc]) delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]+(((((8 delta[FundCol,FCa,FCb]) delta[FundCol,FCc,FCd]) delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi])
+1/4 lambdaSmPmadj[p1,p2,p3,p4] ((((((-2 delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]+(((((2 delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma5[Si,Sl]) gamma5[Sk,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]+(((((2 delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]-(((((2 delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma5[Si,Sj]) gamma5[Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]+(((((8 delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]-(((((8 gamma5[Si,Sl]) gamma5[Sk,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]-(((((8 delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi]+(((((8 gamma5[Si,Sj]) gamma5[Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi])
+1/4 lambdaSpPmadj[p1,p2,p3,p4] ((((((-2 delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]-(((((2 delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma5[Si,Sl]) gamma5[Sk,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]+(((((2 delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]+(((((2 delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma5[Si,Sj]) gamma5[Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]+(((((8 delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]+(((((8 gamma5[Si,Sl]) gamma5[Sk,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]-(((((8 delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi]-(((((8 gamma5[Si,Sj]) gamma5[Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi])
+1/4 lambdaSpPpadj[p1,p2,p3,p4] ((((((-2 delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]-(((((2 delta[FundFlav,Ff,Fi]) delta[FundFlav,Fg,Fh]) gamma5[Si,Sl]) gamma5[Sk,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]+(((((2 delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]+(((((2 delta[FundFlav,Ff,Fg]) delta[FundFlav,Fh,Fi]) gamma5[Si,Sj]) gamma5[Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]-(((((8 delta[Spinor,Si,Sl]) delta[Spinor,Sj,Sk]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]-(((((8 gamma5[Si,Sl]) gamma5[Sk,Sj]) SUNT[Ca,FCa,FCd]) SUNT[Ca,FCc,FCb]) SUNTF[SQCDMFa,Ff,Fi]) SUNTF[SQCDMFa,Fh,Fg]+(((((8 delta[Spinor,Si,Sj]) delta[Spinor,Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi]+(((((8 gamma5[Si,Sj]) gamma5[Sk,Sl]) SUNT[Ca,FCa,FCb]) SUNT[Ca,FCc,FCd]) SUNTF[SQCDMFa,Ff,Fg]) SUNTF[SQCDMFa,Fh,Fi])
)/.{SQCDmu->Unique[SQCDmu],SQCDSn->Unique[SQCDSn],SQCDSm->Unique[SQCDSm],Ca->Unique[Ca],SQCDMFa->Unique[SQCDMFa]}


(* ::Section:: *)
(*Replacement*)


closeChannel={
(*lambdaqqa1[__]->0,*)
lambdaqqa2[__]->0,
lambdaqqa3[__]->0,
(*lambdaqqa4[__]->0,*)
lambdaqqa5[__]->0,
lambdaqqa6[__]->0,
(*lambdaqqa7[__]->0,*)
lambdaqqa8[__]->0
};


closeChannelFourQuark={
lambdaVmA[__]->0,
lambdaVpA[__]->0,
lambdaVmAadj[__]->0,
(*lambdaSigma[__]->0,*)
(*lambdaPion[__]->0,*)
lambdaEta[__]->0,
lambdaA[__]->0,
lambdaSmPmadj[__]->0,
lambdaSpPmadj[__]->0,
lambdaSpPpadj[__]->0
};


rep2point={mq[_]->m,Zq[_]->1,dtZ[_]->0,Rq[_,ss]->Rqs,dRq[_,ss]->dRqs,Rq[_,sr]->Rqr};
rep3point={lambdaqqa1[x___]->lambdaqqa1,
lambdaqqa2[x___]->lambdaqqa2,
lambdaqqa3[x___]->lambdaqqa3,
lambdaqqa4[x___]->lambdaqqa4,
lambdaqqa5[x___]->lambdaqqa5,
lambdaqqa6[x___]->lambdaqqa6,
lambdaqqa7[x___]->lambdaqqa7,
lambdaqqa8[x___]->lambdaqqa8};
rep4point={lambdaVmA[x___]->lambdaVmA,
lambdaVpA[x___]->lambdaVpA,
lambdaVmAadj[x___]->lambdaVmAadj,
lambdaSigma[x___]->lambdaSigma,
lambdaPion[x___]->lambdaPion,
lambdaEta[x___]->lambdaEta,
lambdaA[x___]->lambdaA,
lambdaSmPmadj[x___]->lambdaSmPmadj,
lambdaSpPmadj[x___]->lambdaSpPmadj,
lambdaSpPpadj[x___]->lambdaSpPpadj};


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


(* ::Section:: *)
(*Four  Quark  Vertex*)


Print["Four  Quark  Vertex"];


(* ::Subsection:: *)
(*projector*)


FFTensorVmA[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaVmA[p1,p2,p3,p4]]//Simplify
FFTensorVpA[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaVpA[p1,p2,p3,p4]]//Simplify
FFTensorVmAadj[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaVmAadj[p1,p2,p3,p4]]//Simplify
FFTensorSigma[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaSigma[p1,p2,p3,p4]]//Simplify
FFTensorPion[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaPion[p1,p2,p3,p4]]//Simplify
FFTensorEta[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaEta[p1,p2,p3,p4]]//Simplify
FFTensorA[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaA[p1,p2,p3,p4]]//Simplify
FFTensorSmPmadj[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaSmPmadj[p1,p2,p3,p4]]//Simplify
FFTensorSpPmadj[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaSpPmadj[p1,p2,p3,p4]]//Simplify
FFTensorSpPpadj[FCa_,Si_,Ff_,FCb_,Sj_,Fg_,FCc_,Sk_,Fh_,FCd_,Sl_,Fi_]:=Coefficient[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4],lambdaSpPpadj[p1,p2,p3,p4]]//Simplify


(* ::Subsection:: *)
(*Orthogonality*)


NJL4Fermitracetest[expr_]:=Module[{s},s=(ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk], Sl->Unique[Sl],FCa->Unique[FCa], FCb->Unique[FCb], FCc->Unique[FCc],FCd->Unique[FCd], Ff->Unique[Ff], Fg->Unique[Fg],Fh->Unique[Fh], Fi->Unique[Fi]}];


NJL4FermiFaktortest=FormTrace[NJL4Fermitracetest[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point];


Mtest={Coefficient[NJL4FermiFaktortest,lambdaVmA],Coefficient[NJL4FermiFaktortest,lambdaVpA],Coefficient[NJL4FermiFaktortest,lambdaVmAadj],Coefficient[NJL4FermiFaktortest,lambdaSigma],Coefficient[NJL4FermiFaktortest,lambdaPion],Coefficient[NJL4FermiFaktortest,lambdaEta],Coefficient[NJL4FermiFaktortest,lambdaA],Coefficient[NJL4FermiFaktortest,lambdaSmPmadj],Coefficient[NJL4FermiFaktortest,lambdaSpPmadj],Coefficient[NJL4FermiFaktortest,lambdaSpPpadj]}//Simplify;


FFC1=Solve[Mtest=={1,0,0,0,0,0,0,0,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC2=Solve[Mtest=={0,1,0,0,0,0,0,0,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC3=Solve[Mtest=={0,0,1,0,0,0,0,0,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC4=Solve[Mtest=={0,0,0,1,0,0,0,0,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC5=Solve[Mtest=={0,0,0,0,1,0,0,0,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC6=Solve[Mtest=={0,0,0,0,0,1,0,0,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC7=Solve[Mtest=={0,0,0,0,0,0,1,0,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC8=Solve[Mtest=={0,0,0,0,0,0,0,1,0,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC9=Solve[Mtest=={0,0,0,0,0,0,0,0,1,0},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;
FFC10=Solve[Mtest=={0,0,0,0,0,0,0,0,0,1},{ZVmA,ZVpA,ZVmAadj,ZSigma,ZPion,ZEta,ZA,ZSmPmadj,ZSpPmadj,ZSpPpadj}]//Simplify;


NJL4FermitraceVmA[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC1[[1]])//Simplify)expr ;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceVpA[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC2[[1]])//Simplify)expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceVmAadj[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC3[[1]])//Simplify)expr ;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceSigma[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC4[[1]])//Simplify)expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitracePion[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC5[[1]])//Simplify)expr ;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceEta[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC6[[1]])//Simplify)expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceA[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC7[[1]])//Simplify)expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceSmPmadj[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC8[[1]])//Simplify)expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceSpPmadj[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC9[[1]])//Simplify)expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];
NJL4FermitraceSpPpadj[expr_]:=Module[{s},s=(((ZVmAadj FFTensorVmAadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVmA FFTensorVmA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSigma FFTensorSigma[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZPion FFTensorPion[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZVpA FFTensorVpA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZEta FFTensorEta[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZA FFTensorA[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPmadj FFTensorSpPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSmPmadj FFTensorSmPmadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff]+ZSpPpadj FFTensorSpPpadj[FCd,Sl,Fi,FCc,Sk,Fh,FCb,Sj,Fg,FCa,Si,Ff])/.FFC10[[1]])//Simplify)expr;
s/.{Si->Unique[Si], Sj->Unique[Sj], Sk->Unique[Sk],Sl->Unique[Sl], FCa->Unique[FCa], FCb->Unique[FCb],FCc->Unique[FCc], FCd->Unique[FCd], Ff->Unique[Ff],Fg->Unique[Fg], Fh->Unique[Fh], Fi->Unique[Fi],Ca->Unique[Ca]}];


Print[Simplify[FormTrace[NJL4FermitraceVmA[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaVmA]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceVpA[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaVpA]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceVmAadj[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaVmAadj]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceSigma[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaSigma]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitracePion[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaPion]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceEta[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaEta]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceA[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaA]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceSmPmadj[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaSmPmadj]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceSpPmadj[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaSpPmadj]/.{NRf->2,NRc->3}];
Print[Simplify[FormTrace[NJL4FermitraceSpPpadj[Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCc,Fh,Sk,FCd,Fi,Sl,p1,p2,p3,p4]]/.rep4point]/lambdaSpPpadj]/.{NRf->2,NRc->3}];


(* ::Subsection:: *)
(*diagrams*)


Lambda4qa=(-1)*
(-1)Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCa1,Ff1,Si1,FCb1,Fg1,Sj1,0,0,q,q]Gqfull[FCb1,FCd1,Fg1,Fi1,Sj1,Sl1,q,ss]*
(-1)Gammaq4[FCc,Fh,Sk,FCd,Fi,Sl,FCd1,Fi1,Sl1,FCc1,Fh1,Sk1,0,0,q,q]Gqfull[FCc1,FCa1,Fh1,Ff1,Sk1,Si1,q,sr];


dLambda4qa=(-1)*
(-1)Gammaq4[FCa,Ff,Si,FCb,Fg,Sj,FCa1,Ff1,Si1,FCb1,Fg1,Sj1,0,0,q,q]Gqfull[FCb1,FCd1,Fg1,Fi1,Sj1,Sl1,q,ss]*
(-1)Gammaq4[FCc,Fh,Sk,FCd,Fi,Sl,FCd1,Fi1,Sl1,FCc1,Fh1,Sk1,0,0,q,q]dGqfull[FCc1,FCa1,Fh1,Ff1,Sk1,Si1,q,sr];


Lambda4qb=(1)*
(-1)Gammaq4[FCc,Fh,Sk,FCb,Fg,Sj,FCc1,Fh1,Sk1,FCb1,Fg1,Sj1,0,0,q,q]Gqfull[FCb1,FCd1,Fg1,Fi1,Sj1,Sl1,q,ss]*
(-1)Gammaq4[FCa,Ff,Si,FCd,Fi,Sl,FCd1,Fi1,Sl1,FCa1,Ff1,Si1,0,0,q,q]Gqfull[FCa1,FCc1,Ff1,Fh1,Si1,Sk1,q,sr];


dLambda4qb=(1)*
(-1)Gammaq4[FCc,Fh,Sk,FCb,Fg,Sj,FCc1,Fh1,Sk1,FCb1,Fg1,Sj1,0,0,q,q]Gqfull[FCb1,FCd1,Fg1,Fi1,Sj1,Sl1,q,ss]*
(-1)Gammaq4[FCa,Ff,Si,FCd,Fi,Sl,FCd1,Fi1,Sl1,FCa1,Ff1,Si1,0,0,q,q]dGqfull[FCa1,FCc1,Ff1,Fh1,Si1,Sk1,q,sr];


Lambda4qc=(1/2)*
(-1)Gammaq4[FCa,Ff,Si,FCa1,Ff1,Si1,FCc,Fh,Sk,FCc1,Fh1,Sk1,0,q,0,-q]*Gqfull[FCa1,FCb1,Ff1,Fg1,Si1,Sj1,q,ss]*
Gqfull[FCc1,FCd1,Fh1,Fi1,Sk1,Sl1,-q,sr]*(-1)Gammaq4[FCb1,Fg1,Sj1,FCb,Fg,Sj,FCd1,Fi1,Sl1,FCd,Fi,Sl,q,0,-q,0];


dLambda4qc=(1/2)*
(-1)Gammaq4[FCa,Ff,Si,FCa1,Ff1,Si1,FCc,Fh,Sk,FCc1,Fh1,Sk1,0,q,0,-q]*dGqfull[FCa1,FCb1,Ff1,Fg1,Si1,Sj1,q,ss]*
Gqfull[FCc1,FCd1,Fh1,Fi1,Sk1,Sl1,-q,sr]*(-1)Gammaq4[FCb1,Fg1,Sj1,FCb,Fg,Sj,FCd1,Fi1,Sl1,FCd,Fi,Sl,q,0,-q,0];


dLambda4qa=dLambda4qa/.{FCa1 -> Unique[FCa1],Ff1 -> Unique[Ff1],Si1 -> Unique[Si1],FCb1 -> Unique[FCb1],
Fg1 -> Unique[Fg1], Sj1 -> Unique[Sj1], FCd1 -> Unique[FCd1], Fi1 -> Unique[Fi1], Sl1 -> Unique[Sl1],
FCc1 -> Unique[FCc1], Fh1 -> Unique[Fh1], Sk1 -> Unique[Sk1]};
dLambda4qb=dLambda4qb/.{FCa1 -> Unique[FCa1],Ff1 -> Unique[Ff1],Si1 -> Unique[Si1],FCb1 -> Unique[FCb1],
Fg1 -> Unique[Fg1], Sj1 -> Unique[Sj1], FCd1 -> Unique[FCd1], Fi1 -> Unique[Fi1], Sl1 -> Unique[Sl1],
FCc1 -> Unique[FCc1], Fh1 -> Unique[Fh1], Sk1 -> Unique[Sk1]};
dLambda4qc=dLambda4qc/.{FCa1 -> Unique[FCa1],Ff1 -> Unique[Ff1],Si1 -> Unique[Si1],FCb1 -> Unique[FCb1],
Fg1 -> Unique[Fg1], Sj1 -> Unique[Sj1], FCd1 -> Unique[FCd1], Fi1 -> Unique[Fi1], Sl1 -> Unique[Sl1],
FCc1 -> Unique[FCc1], Fh1 -> Unique[Fh1], Sk1 -> Unique[Sk1]};
(*lambdaAlgAE=dLambda4qa+dLambda4qb+dLambda4qc;*)


NJL4FermiVmAa=NJL4FermitraceVmA[dLambda4qa/.rep4point];
NJL4FermiVpAa=NJL4FermitraceVpA[dLambda4qa/.rep4point];
NJL4FermiVmAadja=NJL4FermitraceVmAadj[dLambda4qa/.rep4point];
NJL4FermiSigmaa=NJL4FermitraceSigma[dLambda4qa/.rep4point];
NJL4FermiPiona=NJL4FermitracePion[dLambda4qa/.rep4point];
NJL4FermiEtaa=NJL4FermitraceEta[dLambda4qa/.rep4point];
NJL4FermiAa=NJL4FermitraceA[dLambda4qa/.rep4point];
NJL4FermiSmPmadja=NJL4FermitraceSmPmadj[dLambda4qa/.rep4point];
NJL4FermiSpPmadja=NJL4FermitraceSpPmadj[dLambda4qa/.rep4point];
NJL4FermiSpPpadja=NJL4FermitraceSpPpadj[dLambda4qa/.rep4point];

NJL4FermiVmAb=NJL4FermitraceVmA[dLambda4qb/.rep4point];
NJL4FermiVpAb=NJL4FermitraceVpA[dLambda4qb/.rep4point];
NJL4FermiVmAadjb=NJL4FermitraceVmAadj[dLambda4qb/.rep4point];
NJL4FermiSigmab=NJL4FermitraceSigma[dLambda4qb/.rep4point];
NJL4FermiPionb=NJL4FermitracePion[dLambda4qb/.rep4point];
NJL4FermiEtab=NJL4FermitraceEta[dLambda4qb/.rep4point];
NJL4FermiAb=NJL4FermitraceA[dLambda4qb/.rep4point];
NJL4FermiSmPmadjb=NJL4FermitraceSmPmadj[dLambda4qb/.rep4point];
NJL4FermiSpPmadjb=NJL4FermitraceSpPmadj[dLambda4qb/.rep4point];
NJL4FermiSpPpadjb=NJL4FermitraceSpPpadj[dLambda4qb/.rep4point];

NJL4FermiVmAc=NJL4FermitraceVmA[dLambda4qc/.rep4point];
NJL4FermiVpAc=NJL4FermitraceVpA[dLambda4qc/.rep4point];
NJL4FermiVmAadjc=NJL4FermitraceVmAadj[dLambda4qc/.rep4point];
NJL4FermiSigmac=NJL4FermitraceSigma[dLambda4qc/.rep4point];
NJL4FermiPionc=NJL4FermitracePion[dLambda4qc/.rep4point];
NJL4FermiEtac=NJL4FermitraceEta[dLambda4qc/.rep4point];
NJL4FermiAc=NJL4FermitraceA[dLambda4qc/.rep4point];
NJL4FermiSmPmadjc=NJL4FermitraceSmPmadj[dLambda4qc/.rep4point];
NJL4FermiSpPmadjc=NJL4FermitraceSpPmadj[dLambda4qc/.rep4point];
NJL4FermiSpPpadjc=NJL4FermitraceSpPpadj[dLambda4qc/.rep4point];


flow4={
NJL4FermiVmAa,NJL4FermiVpAa,NJL4FermiVmAadja,NJL4FermiSigmaa,NJL4FermiPiona,NJL4FermiEtaa,NJL4FermiAa,NJL4FermiSmPmadja,NJL4FermiSpPmadja,NJL4FermiSpPpadja,
NJL4FermiVmAb,NJL4FermiVpAb,NJL4FermiVmAadjb,NJL4FermiSigmab,NJL4FermiPionb,NJL4FermiEtab,NJL4FermiAb,NJL4FermiSmPmadjb,NJL4FermiSpPmadjb,NJL4FermiSpPpadjb,
NJL4FermiVmAc,NJL4FermiVpAc,NJL4FermiVmAadjc,NJL4FermiSigmac,NJL4FermiPionc,NJL4FermiEtac,NJL4FermiAc,NJL4FermiSmPmadjc,NJL4FermiSpPmadjc,NJL4FermiSpPpadjc
};


Print["FormTracer begin"]
LaunchKernels[];
DistributeDefinitions[FormTrace];
flow4form=ParallelMap[(FormTrace[#])&,flow4];
Print["FormTracer OK"]
Export["flow4formnop.wdx",flow4form,"WDX"];
