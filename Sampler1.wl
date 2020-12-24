(* ::Package:: *)

BeginPackage["Sampler1`"]
hmc::usage = "qs=hmc[U,dU,ddU,Dim,BURNIN,EPISODE,vanilla,switch]";
HessianH::usage="HessianH[f,x]";
GradientG::usage="GradientG[f,x]";
psrf::usage="psrf[QS]";
CHAINS::usage = "mcmc chains, or particles";
dt10::usage = "time step for hmc";
dt20::usage = "time step for rhmc";
STEPS::usage="steps"
INTERVAL::usage="INTERVAL"


Begin["`Private`"]

CHAINS=10;
STEPS=10;
dt10=0.00001;
dt20=0.01;
INTERVAL=1001;
hmc[U_,dU_,ddU_,Dim_,BURNIN_,EPISODE_,vanilla0_,switch_]:=Module[{qAll,pAll,Utotal,Ktotal,Htotal,s,S,AS,ES,KtotalNew,dt,p,q0,UE,\[Alpha],q,j,Htotal1,Htotal2,
ND=NormalDistribution[0,1],
UD=UniformDistribution[],
QS=List[],vanilla=vanilla0,dt1=dt10,dt2=dt20,ACS={}},
pAll=RandomVariate[ND,{CHAINS,Dim}];
qAll=RandomVariate[ND,{CHAINS,Dim}];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotal1=2Utotal;
Htotal2=2Utotal;
For[j=1,j\[LessSlantEqual]EPISODE,j++,
pAll=RandomVariate[ND,{CHAINS,Dim}];
KtotalNew = Sum[If[vanilla,pAll[[i]].pAll[[i]],pAll[[i]]. LinearSolve[Apply[ddU,qAll[[i]]],pAll[[i]]]]/2,{i,1,CHAINS}];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotal=If[vanilla,Htotal1,Htotal2];
dt=If[vanilla,dt1,dt2];
Ktotal=Htotal-Utotal;
pAll=pAll Sqrt[Abs[Ktotal/KtotalNew]];
ES={};
AS={};
For[i=1,i<=CHAINS,i++,
p=pAll[[i]];
q=qAll[[i]];
UE={Apply[ U,q]};
q0=q;
For[s=1;,s<=STEPS,s++,
p=p-dt Apply[dU,q];
q=q+dt If[vanilla,p,LinearSolve[Apply[ddU,q],p]];
UE=Append[UE,Apply[ U,q]]];
ES=Append[ES,UE];
\[Alpha]=Exp[Clip[Apply[U,q0]-Apply[U,q],{-20,0}]];
AS=Append[AS,N[\[Alpha]]];
If[\[Alpha]<RandomVariate[UD],q=q0];
qAll[[i]]=q;
If[j>BURNIN,QS=Append[QS,q]]];
s=Union[Flatten[Table[Ordering[ES[[i]],1],{i,1,CHAINS}]]];
S =Union[Flatten[Table[Ordering[ES[[i]],-1],{i,1,CHAINS}]]];
If[True  && Mod[j,INTERVAL]==0,Print[j," ",Ktotal," ",KtotalNew," ",Mean[AS]," ",StandardDeviation[AS]," ",Utotal," ", Htotal1," ",Htotal2," ",dt1," ",dt2," ", vanilla," ",s," ",S]];
If[j<BURNIN,
If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt 1.1];
If[s=={1} && S=={STEPS+1},dt=dt/1.1];
If[Ktotal>0&&KtotalNew>0,
If[Utotal<0&&Htotal>Utotal+1000,Htotal=Utotal+1000,If[Utotal>0&&Htotal>100Utotal,Htotal=100Utotal,
If[Mean[AS]>.9,Ktotal=Ktotal 1.1,
If[Mean[AS]<.1,Ktotal=Ktotal/ 1.1]];Htotal=Ktotal+Utotal]]];
If[vanilla,
Htotal1=Htotal;
dt1=dt,
Htotal2=Htotal;
dt2=dt]];
If[switch,vanilla=Not[vanilla]]];
QS];


psrf[QS1_]:=Module[{QQ,\[Theta],\[Sigma]2,B,W,V,Rs},
QQ=Table[Table[QS1[[i]],{i,j,50000+j-1,10}],{j,1,10}];
\[Theta]=Table[Mean[QQ[[j]]],{j,1,10}];
\[Sigma]2=Table[Variance[QQ[[j]]],{j,1,10}];
B=(Total[Table[\[Theta][[j]]-Mean[\[Theta]],{j,1,10}]^2] 5000)/(10-1);
W=Mean[\[Sigma]2];
V=(5000-1)/5000 W+(10+1)/(5000 10) B;
Rs=Sqrt[V/W];
Rs]

HessianH[f_,x_List?VectorQ]:=N[D[f,{x,2}]];
GradientG[f_,x_List?VectorQ]:=N[D[f,{x,1}]];


End[]


EndPackage[]
