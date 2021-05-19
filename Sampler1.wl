(* ::Package:: *)

BeginPackage["Sampler1`"]
hmc::usage = "qs=hmc[U,dU,ddU,Dim,BURNIN,EPISODE,vanilla,switch,qinit]";
hmc1::usage = "qs=hmc[U,dU,ddU,Dim,BURNIN,EPISODE,vanilla,switch,qinit]";
hmc0::usage = "qs=hmc[U,dU,ddU,Dim,BURNIN,EPISODE,vanilla]";
HessianH::usage="HessianH[f,x]";
GradientG::usage="GradientG[f,x]";
psrf::usage="psrf[QS]";
CHAINS::usage = "mcmc chains, or particles";
dt10::usage = "time step for hmc";
dt20::usage = "time step for hmc";
STEPS::usage="steps"
INTERVAL::usage="INTERVAL"
outbnd::usage="outbnd"
DECAYDT::usage="DECAYDT";
DECAYENERGY::usage="DECAYENERGY";
KU1::usage="KU1";
KU2::usage="KU2";
KU::usage="KU";
LOWAP::usage="LOWAP";
HIGHAP::usage="HIGHAP";
switch0::usage="switch0";
diag::usage="diag";
forcetune::usage="forcetune";
mass::usage="mass";
callback::usage="callback";

Begin["`Private`"]
CHAINS=3;
STEPS=3;
dt10=0.00000000000001;
dt20=0.00000000000001;
INTERVAL=1001;
outbnd[q_]:=False;
callback[j_]:=None;
DECAYDT=.1;
DECAYENERGY=.1;
KU=1000;
LOWAP=0.1;
HIGHAP=0.9;
switch0=100;
diag=False;

hmc2[U_,dU_,ddU_,Dim_,BURNIN_,EPISODE_,vanilla0_,switch_,qinit_]:=Module[{qAll,pAll,Utotal,Ktotal,hi,lo,decayenergy=DECAYENERGY,decaydt=DECAYDT,decayfactor=100^(-1/BURNIN),dt1=dt10,dt2=dt20,Htotal,UUtotal,s,S,s1,AS,m,KtotalNew,dt,p,q0,\[Alpha],q,j,bad,Htotal1,Htotal2,i,q1,ND=NormalDistribution[0,1],
UD=UniformDistribution[],QS=List[],vanilla=vanilla0,KtotalNews,ACS={}},
If[qinit!={},qAll=qinit,qAll=RandomVariate[ND,{CHAINS,Dim}]];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
UUtotal=Utotal;
Htotal1=If[Utotal>0,2Utotal,100];
Htotal2=If[Utotal>0,2Utotal,100];
For[j=1,j\[LessSlantEqual]EPISODE,j++,
pAll=RandomVariate[ND,{CHAINS,Dim}];
KtotalNews = Table[If[vanilla,pAll[[i]].pAll[[i]],pAll[[i]].If[diag,pAll[[i]]/Apply[ddU,qAll[[i]]],LinearSolve[Apply[ddU,qAll[[i]]],pAll[[i]]]]]/2,{i,1,CHAINS}];
KtotalNew = Total[KtotalNews];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotal=If[vanilla,Htotal1,Htotal2];dt=If[vanilla,dt1,dt2];Ktotal=Htotal-Utotal;pAll=pAll Sqrt[Abs[Ktotal/KtotalNew]];AS={};
For[i=1,i<=CHAINS,i++,bad=False;p=pAll[[i]];q=qAll[[i]];q0=q;For[s=1;,s<=STEPS,s++,p=p-dt Apply[dU,q];q1=q;q=q+dt If[vanilla,p,If[diag,p/Apply[ddU,q],LinearSolve[Apply[ddU,q],p]]];If[outbnd[q],bad=True;Break[]]];\[Alpha]=If[bad,0,Exp[Clip[Apply[U,q0]-Apply[U,q],{-20,0}]]];AS=Append[AS,N[\[Alpha]]];If[bad || \[Alpha]<RandomVariate[UD],q=q0];qAll[[i]]=q;If[j>BURNIN,QS=Append[QS,q]]];If[Mod[j,INTERVAL]==0,Print[Row[{j,Ktotal,KtotalNew,Mean[AS],StandardDeviation[AS],Utotal,Htotal1,Htotal2,dt1,dt2,vanilla},"     "]]];
If[j<BURNIN,
If[Mean[AS]>HIGHAP,dt=dt (1+decaydt)];
If[Mean[AS]<LOWAP,dt=dt/(1+decaydt)];
inc=Mean[AS]>HIGHAP;
dec=Mean[AS]<LOWAP;
If[Htotal<Utotal,If[inc,dec=True;inc=False];If[dec,inc=True;dec=False]];
If[inc,If[Htotal>0,Htotal=Htotal (1+decayenergy),Htotal=Htotal/(1+decayenergy)]];
If[dec,If[Htotal>0,Htotal=Htotal/(1+decayenergy),Htotal=Htotal (1+decayenergy)]];
If[vanilla,Htotal1=Htotal;dt1=dt,Htotal2=Htotal;dt2=dt]];
If[switch,vanilla=Not[vanilla]]];QS];

hmc[U_,dU_,ddU_,Dim_,BURNIN_,EPISODE_,vanilla0_,switch_,qinit_]:=Module[{qAll,pAll,hess,UE,ES,Utotal,Ktotal,hi,lo,decayenergy=DECAYENERGY,decaydt=DECAYDT,dt1=dt10,decayfactor=10^(-1/BURNIN),dt2=dt20,Htotal,ii=1,s,S,AS,KtotalNew,dt,p,q0,\[Alpha],q,j,bad,Htotal1,Htotal2,i,q1,ND=NormalDistribution[0,1],UD=UniformDistribution[],QS=ConstantArray[0,{CHAINS( EPISODE-BURNIN),Dim}],vanilla=vanilla0,ACS={}},
If[qinit!={},qAll=qinit,qAll=RandomVariate[ND,{CHAINS,Dim}]];Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];Htotal1=2Utotal;Htotal2=2Utotal;
For[j=1,j\[LessSlantEqual]EPISODE,j++,
pAll=RandomVariate[ND,{CHAINS,Dim}];KtotalNew = Sum[If[vanilla,pAll[[i]].pAll[[i]],pAll[[i]].If[diag,pAll[[i]]/Apply[ddU,qAll[[i]]],LinearSolve[Apply[ddU,qAll[[i]]],pAll[[i]],Method->"Krylov"]]]/2,{i,1,CHAINS}];Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];Htotal=If[vanilla,Htotal1,Htotal2];dt=If[vanilla,dt1,dt2];Ktotal=Htotal-Utotal;pAll=pAll Sqrt[Abs[Ktotal/KtotalNew]];AS={};ES={};For[i=1,i<=CHAINS,i++,
bad=False;
p=pAll[[i]];
q=qAll[[i]];
q0=q;
If[j<=BURNIN,UE={Apply[U,q]}];Do[p=p-dt Apply[dU,q];q1=q;q=q+dt If[vanilla,p,If[diag,p/Apply[ddU,q],LinearSolve[Apply[ddU,q],p,Method->"Krylov"]]];
If[j<=BURNIN,UE=Append[UE,Apply[U,q]]],STEPS];If[outbnd[q],bad=True];
If[j<=BURNIN,ES=Append[ES,UE]];\[Alpha]=If[bad,0,Exp[Clip[Apply[U,q0]-Apply[U,q],{-20,0}]]];AS=Append[AS,N[\[Alpha]]];If[\[Alpha]<RandomVariate[UD],q=q0];qAll[[i]]=q;If[j>BURNIN,QS[[ii,;;]]=q;ii=ii+1]];If[Mod[j,INTERVAL]==0,Print[Row[{j,Ktotal,KtotalNew,Mean[AS],StandardDeviation[AS],Utotal,Htotal1,Htotal2,dt1,dt2,vanilla,s,S},"     "]]];If[j<=BURNIN,
s=Union[Flatten[Table[Ordering[ES[[i]],1],{i,1,CHAINS}]]];S=Union[Flatten[Table[Ordering[ES[[i]],-1],{i,1,CHAINS}]]];If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt (1+decaydt)];If[s=={1}&&S=={STEPS+1},dt=dt/(1+decaydt)];hi=Mean[AS]>HIGHAP;lo=Mean[AS]<LOWAP;If[hi,Htotal=(Htotal-Utotal)(1+decayenergy)+Utotal,If[lo,Htotal=(Htotal-Utotal)/(1+decayenergy)+Utotal]]];If[vanilla,
Htotal1=Htotal;dt1=dt,
Htotal2=Htotal;dt2=dt];
If[switch && j>switch0,
vanilla=Not[vanilla]];callback[j]];QS];

hmc1[U_,dU_,ddU_,Dim_,BURNIN_,EPISODE_,vanilla0_,switch_,qinit_]:=Module[{qAll,pAll,decaydt=DECAYDT,decayenergy=DECAYENERGY,decayfactor=1000^(-1/BURNIN),updt,Utotal,Ktotal,Htotal,s,S,AS,ES,KtotalNew,dt,p,q0,UE,\[Alpha],q,j,bad,Htotal1,Htotal2,i,q1,ND=NormalDistribution[0,1],UD=UniformDistribution[],QS=List[],vanilla=vanilla0,dt1=dt10,dt2=dt20,KtotalNews,ACS={}},If[qinit!={},qAll=qinit,qAll=RandomVariate[ND,{CHAINS,Dim}]];
Utotal=Sum[Apply[U,qAll[[i]]],{i,1,CHAINS}];
Htotal1=2Utotal;
Htotal2=2Utotal;
For[j=1,j\[LessSlantEqual]EPISODE,j++,
(*If[j<100,Htotal1=2Abs[Utotal];Htotal2=2Abs[Utotal]];*)
pAll=RandomVariate[ND,{CHAINS,Dim}];
KtotalNew = Sum[ If[vanilla,pAll[[i]].pAll[[i]],pAll[[i]].If[diag,pAll[[i]]/Apply[ddU,qAll[[i]]],LinearSolve[Apply[ddU,qAll[[i]]],pAll[[i]]]]]/2,{i,1,CHAINS}];
Utotal=Sum[Apply[U,qAll[[i]]],{i,1,CHAINS}];Htotal=If[vanilla,Htotal1,Htotal1];dt=If[vanilla,dt1,dt2];Ktotal=Htotal-Utotal;pAll=pAll Sqrt[Abs[Ktotal/KtotalNew]];ES={};AS={};For[i=1,i<=CHAINS,i++,
p=pAll[[i]];
q=qAll[[i]];
UE={Apply[U,q]};
q0=q;
bad=False;
For[s=1;,s<=STEPS,s++,
If[Not[bad],
p=p-dt Apply[dU,q];
q1=q;
q=q+dt If[vanilla,p,LinearSolve[Apply[ddU,q],p]];
If[outbnd[q],q=q1;bad=True]];
UE=Append[UE,Apply[U,q]]];
ES=Append[ES,UE];
\[Alpha]=If[bad,0,Exp[Clip[Apply[U,q0]-Apply[U,q],{-20,0}]]];
AS=Append[AS,N[\[Alpha]]];
If[\[Alpha]<RandomVariate[UD],q=q0];
qAll[[i]]=q;
If[j>BURNIN,QS=Append[QS,q]]];
s=Union[Flatten[Table[Ordering[ES[[i]],1],{i,1,CHAINS}]]];
S=Union[Flatten[Table[Ordering[ES[[i]],-1],{i,1,CHAINS}]]];
If[Mod[j,INTERVAL]==0,Print[j," ",Ktotal," ",KtotalNew," ",Mean[AS]," ",StandardDeviation[AS]," ",Utotal," ",Htotal1," ",Htotal2," ",dt1," ",dt2," ",vanilla," ",s," ",S," ",decaydt]];

(*If[RandomChoice[{True,False}],*)

If[j<BURNIN ,
If[Htotal>0,
If[Mean[AS]>HIGHAP,Htotal=Htotal (1+decayenergy)];
If[Mean[AS]<LOWAP,Htotal=Htotal/(1+decayenergy)]];
If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt (1+decaydt)];
If[s=={1}&&S=={STEPS+1},dt=dt/(1+decaydt)];
decaydt=decaydt decayfactor;
(*decayenergy=decayenergy decayfactor;*)
If[vanilla,dt1=dt,dt2=dt];
If[vanilla,Htotal1=Htotal,Htotal1=Htotal]];
(*Htotal=Utotal+Ktotal;*)
(*If[Htotal<KU2 Utotal,Htotal=KU2 Utotal,If[Htotal>KU1 Utotal,Htotal=KU1 Utotal]];*)
If[switch,vanilla=Not[vanilla]]];QS];

HessianH[f_,x_List?VectorQ]:=N[D[f,{x,2}]];
GradientG[f_,x_List?VectorQ]:=N[D[f,{x,1}]];


hmc0[U_,dU_,ddU_,Dim_,BURNIN_,EPISODE_,vanilla_,qinit_]:=Module[{qAll,pAll,Utotal,Ktotal,s,S,AS,ES,p,p0,q0,UE,\[Alpha],q,j,bad,q1,i,ND=NormalDistribution[0,1],UD=UniformDistribution[],QS=List[],dt=dt10,ACS={}},
If[qinit!={},qAll=qinit,qAll=RandomVariate[ND,{CHAINS,Dim}]];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];For[j=1,j\[LessSlantEqual]EPISODE,j++,
pAll=RandomVariate[ND,{CHAINS,Dim}];Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];ES={};AS={};For[i=1,i<=CHAINS,i++,bad=False;
p=pAll[[i]];
q=qAll[[i]];
UE={Apply[ U,q]};
q0=q;
p0=p;
For[s=1;,s<=STEPS,s++,
p=p-dt Apply[dU,q];
q1=q;
q=q+dt If[vanilla,p,LinearSolve[Apply[ddU,q],p]];
If[outbnd[q],q=q1;bad=True];
UE=Append[UE,Apply[ U,q]]];
ES=Append[ES,UE];
\[Alpha]=Exp[Clip[Apply[U,q0]-Apply[U,q]+If[vanilla,p0.p0-p.p,p0. LinearSolve[Apply[ddU,q0],p0]-p. LinearSolve[Apply[ddU,q],p]]/2,{-20,0}]];AS=Append[AS,N[\[Alpha]]];If[\[Alpha]<RandomVariate[UD] || bad,q=q0];qAll[[i]]=q;If[j>BURNIN,QS=Append[QS,q]]];s=Union[Flatten[Table[Ordering[ES[[i]],1],{i,1,CHAINS}]]];S =Union[Flatten[Table[Ordering[ES[[i]],-1],{i,1,CHAINS}]]];If[True && Mod[j,INTERVAL]==0,Print[j," ",Mean[AS]," ",StandardDeviation[AS]," ",Utotal," ", dt," ", vanilla," ",s," ",S]];If[j<BURNIN,If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt 1.01];If[s=={1} && S=={STEPS+1},dt=dt/1.01]]];
QS];


psrf[QS_]:=Module[{QQ,\[Theta],\[Sigma]2,B,W,V,Rs,len=Length[QS],chains=CHAINS},
QQ=Table[QS[[i;;Length[QS];;CHAINS,;;]],{i,1,CHAINS}];
\[Theta]=Table[Mean[QQ[[j]]],{j,1,chains}];
\[Sigma]2=Table[Variance[QQ[[j]]],{j,1,chains}];
B=(Total[Table[\[Theta][[j]]-Mean[\[Theta]],{j,1,chains}]^2] len/chains)/(chains-1);
W=Mean[\[Sigma]2];
V=(len/chains-1)/(len/chains) W+(chains+1)/(len/chains chains) B;
Rs=Sqrt[V/W];
Rs]


End[]


EndPackage[]
