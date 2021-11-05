(* ::Package:: *)

BeginPackage["Sampler2`"]
hmc2::usage = "hmc2[U,dU,K,dKdp,Dim,BURNIN,EPISODE,qinit,levels0]";
(*hmc13::usage = "qs=hmc[U,dU,ddU,Dim,BURNIN,EPISODE,vanilla,switch,qinit]";*)
HessianH::usage="HessianH[f,x]";
GradientG::usage="GradientG[f,x]";
psrf::usage="psrf[QS]";
CHAINS::usage = "mcmc chains, or particles";
dt0::usage = "time step for hmc";
(*dt30::usage = "time step for hmc";*)
STEPS::usage="steps"
(*setchain::usage="setchain"*)
INTERVAL::usage="INTERVAL"
outbnd::usage="outbnd"
DECAYDT::usage="DECAYDT";
DECAYENERGY::usage="DECAYENERGY";
LOWLEVEL::usage="LOWAP";
HIGHLEVEL::usage="HIGHAP";

Begin["`Private`"]
CHAINS=3;
STEPS=3;
dt0=10^-6;
INTERVAL=1001;
outbnd[q_]:=False;
DECAYDT=.1;
DECAYENERGY=.1;
LOWLEVEL=0.1;
HIGHLEVEL=0.9;

hmc2[U_,dU_,K_,dKdp_,Dim_,BURNIN_,EPISODE_,qinit_,levels0_]:=Module[{
qAll,pAll,UE,ES,Utotal,Ktotal,hi,lo,Htotal,Htotals,s,S,AS,KtotalNew,dt,p,q0,\[Alpha],q,j,bad,i,q1,decayenergy=DECAYENERGY,decaydt=DECAYDT,ii=1, k0,h0,k1,h1,p1,QS=ConstantArray[0,{CHAINS( EPISODE-BURNIN),Dim}],ACS={},dts=Table[dt0,Length[levels0]],level,es,h,pos1=True,pos2,levels=levels0},
If[qinit!={},
qAll=qinit,
qAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}]];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotals=Table[If[Utotal==0,1,2Utotal],Length[levels]];
level=1;
For[j=1,j\[LessSlantEqual]EPISODE,j++,
(*If[Mod[j,Length[levels0]]==1,pos1=Not[pos1]];*)
(*pos2=pos1;*)
pAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}];
KtotalNew = Sum[Apply[K,Append[Join[pAll[[i]],qAll[[i]]],levels[[level]]]],{i,1,CHAINS}];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotal=Htotals[[level]];
dt=dts[[level]];
Ktotal=Htotal-Utotal;
pAll=pAll Sqrt[Abs[Ktotal/KtotalNew]];
AS={};ES={};
For[i=1,i<=CHAINS,i++,
bad=False;
p=pAll[[i]];
q=qAll[[i]];
q0=q;
If[j<=BURNIN,UE={Apply[U,q]}];
k0=Apply[K,Append[Join[p,q],levels[[level]]]];
h0=Apply[U,q]+k0;
Do[
q1=q;
p1=p;
q=q+dt Apply[dKdp,Append[Join[p,q],levels[[level]]]];
p=p-dt Apply[dU,q];
k1=Apply[K,Append[Join[p,q],levels[[level]]]];
h1=Apply[U,q]+k1;
p=p Sqrt[Abs[(k1-(h1-h0))/k1]];

If[outbnd[q],bad=True;q=q1];
If[j<=BURNIN,UE=Append[UE,Apply[U,q]]],STEPS];
If[j<=BURNIN,ES=Append[ES,UE]];
\[Alpha]=If[bad,0,Exp[Clip[Apply[U,q0]-Apply[U,q],{-10,0}]]];
AS=Append[AS,N[\[Alpha]]];
If[\[Alpha]<RandomVariate[UniformDistribution[]],q=q0];
qAll[[i]]=q;
If[j>BURNIN,
QS[[ii,;;]]=q;
ii=ii+1]];
If[j<=BURNIN,
s=Union[Flatten[Table[Ordering[ES[[i]],1],{i,1,CHAINS}]]];
S=Union[Flatten[Table[Ordering[ES[[i]],-1],{i,1,CHAINS}]]];
If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt (1+decaydt)];
If[s=={1}&&S=={STEPS+1},dt=dt/(1+decaydt)];
hi=Quantile[AS,.5]>HIGHLEVEL;
lo=Quantile[AS,.5]<LOWLEVEL;
If[hi||lo,Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}]];
If[hi,Htotal=(Htotal-Utotal)(1+decayenergy)+Utotal,If[lo,Htotal=(Htotal-Utotal)/(1+decayenergy)+Utotal]];
dts[[level]]=dt;
Htotals[[level]]=Htotal;
];
If[Mod[j,INTERVAL]==0,Print[Row[{j,Ktotal,KtotalNew,Mean[AS],StandardDeviation[AS],Utotal,level,s,S,dts,Htotals,dt},"     "]]];
level=Mod[level,Length[levels]]+1];QS];












HessianH[f_,x_List?VectorQ]:=N[D[f,{x,2}]];
GradientG[f_,x_List?VectorQ]:=N[D[f,{x,1}]];


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



