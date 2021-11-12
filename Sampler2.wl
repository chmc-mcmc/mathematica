(* ::Package:: *)

BeginPackage["Sampler2`"]
hmc2::usage = "hmc2[U,dU,K,dKdp,Dim,BURNIN,EPISODE,qinit,levels0]";
hmc3::usage = "hmc2[U,dU,K,dKdp,dKdq,Dim,BURNIN,EPISODE,qinit,levels0]";
dKdq0::usage="dKdq0[p,q,r0,ddU,dddU]";
get\[CapitalSigma]::usage="get\[CapitalSigma][q,r,ddU]";
HessianH::usage="HessianH[f,x]";
GradientG::usage="GradientG[f,x]";
D3::usage="D3[f,x]";
psrf::usage="psrf[QS]";
CHAINS::usage = "mcmc chains, or particles";
dt0::usage = "time step for hmc";
STEPS::usage="steps"
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
qAll,pAll,UE,ES,Utotal,Ktotal,hi,lo,Htotal,Htotals,s,S,AS,KtotalNew,dt,p,q0,\[Alpha],q,j,bad,switch1,switch=Table[True,Length[levels0]],i,q1,decayenergy=DECAYENERGY,decaydt=DECAYDT,ii=1, k0,h0,k1,h1,p1,QS=ConstantArray[0,{CHAINS( EPISODE-BURNIN),Dim}],ACS={},dts=Table[dt0,Length[levels0]],level,es,h,pos1=True,pos2,levels=levels0},
If[qinit!={},
qAll=qinit,
qAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}]];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotals=Table[If[Utotal==0,1,2Utotal],Length[levels]];
level=1;
For[j=1,j\[LessSlantEqual]EPISODE,j++,
(*switch[[level]]=Not[switch[[level]]];
switch1=switch[[level]];*)
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

(*p=p- dt .5 Apply[dU,q];*)

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

(*q=q+dt Apply[dKdp,Append[Join[p,q],levels[[level]]]];
p=p- dt .5 Apply[dU,q];
*)

(*If[switch1&&levels[[level]]>0,q=2q0-q];*)
(*switch1=Not[switch1];*)
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
If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt(1+decaydt)];
If[s=={1}&&S=={STEPS+1},dt=dt/(1+decaydt)];
hi=Quantile[AS,.5]>HIGHLEVEL;
lo=Quantile[AS,.5]<LOWLEVEL;
If[hi||lo,Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}]];
If[hi,Htotal=(Htotal-Utotal)(1+decayenergy)+Utotal,If[lo,Htotal=(Htotal-Utotal)/(1+decayenergy)+Utotal]];
dts[[level]]=dt;
Htotals[[level]]=Htotal;
];
If[Mod[j,INTERVAL]==0,Print[Row[{j,Ktotal,KtotalNew,Mean[AS],StandardDeviation[AS],Utotal,level,s,S,dts,Htotals,dt,hi,lo},"     "]]];
level=Mod[level,Length[levels]]+1];QS];

hmc3[U_,dU_,K_,dKdp_,dKdq_,Dim_,BURNIN_,EPISODE_,qinit_,levels0_]:=Module[{
qAll,pAll,UE,ES,Utotal,Ktotal,hi,lo,Htotal,Htotals,s,S,AS,KtotalNew,dt,p,q0,\[Alpha],q,j,bad,switch1,k0,h0,k1,h1,i,q1,decayenergy=DECAYENERGY,decaydt=DECAYDT,ii=1, QS=ConstantArray[0,{CHAINS( EPISODE-BURNIN),Dim}],ACS={},dts=Table[dt0,Length[levels0]],level,es,h,pos1=True,pos2,levels=levels0},
If[qinit!={},
qAll=qinit,
qAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}]];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotals=Table[If[Utotal==0,1,2Utotal],Length[levels]];
level=1;
For[j=1,j\[LessSlantEqual]EPISODE,j++,
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

(*k0=Apply[K,Append[Join[p,q],levels[[level]]]];
h0=Apply[U,q]+k0;
*)

p=p- dt .5 (Apply[dU,q]+Apply[dKdq,Append[Join[p,q],levels[[level]]]]);

Do[
q1=q;
q=q+dt Apply[dKdp,Append[Join[p,q],levels[[level]]]];
p=p-dt (Apply[dU,q]+Apply[dKdq,Append[Join[p,q],levels[[level]]]]);

(*k1=Apply[K,Append[Join[p,q],levels[[level]]]];
h1=Apply[U,q]+k1;
p=p Sqrt[Abs[(k1-(h1-h0))/k1]];
*)

If[outbnd[q],bad=True;q=q1];
If[j<=BURNIN,UE=Append[UE,Apply[U,q]]],STEPS];

(*q=q+dt Apply[dKdp,Append[Join[p,q],levels[[level]]]];*)
(*p=p- dt .5 (Apply[dU,q]+Apply[dKdq,Append[Join[p,q],levels[[level]]]]);*)

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
If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt(1+decaydt)];
If[s=={1}&&S=={STEPS+1},dt=dt/(1+decaydt)];
hi=Quantile[AS,.5]>HIGHLEVEL;
lo=Quantile[AS,.5]<LOWLEVEL;
If[hi||lo,Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}]];
If[hi,Htotal=(Htotal-Utotal)(1+decayenergy)+Utotal,If[lo,Htotal=(Htotal-Utotal)/(1+decayenergy)+Utotal]];
dts[[level]]=dt;
Htotals[[level]]=Htotal;
];
If[Mod[j,INTERVAL]==0,Print[Row[{j,Ktotal,KtotalNew,Mean[AS],StandardDeviation[AS],Utotal,level,s,S,dts,Htotals,dt,hi,lo,1},"     "]]];
level=Mod[level,Length[levels]]+1];QS];


(*https://mathoverflow.net/questions/229425/derivative-of-eigenvectors-of-a-matrix-with-respect-to-its-components*)
dKdq0[p_,q_,r0_,ddU_,dddU_]:=Module[{r=-r0,\[Lambda]d,dim,eig,v,e,y,Bdot,vv,dudq},
\[Lambda]d=r Abs[e]^(r-1);
dim=Length[p];
eig=Eigensystem[Apply[ddU,q]];
v=Normalize/@eig[[2]];
e=eig[[1]];
y=v.p;
Bdot=Apply[dddU,q];vv=Table[Table[Sum[If[i!=j,1/(e[[i]]-e[[j]]) (v[[j]].(Bdot[[k]].v[[i]]))v[[j]],0],{j,1,dim}],{i,1,dim}],{k,1,dim}];dudq=Table[p.(Transpose[vv[[i]]].(DiagonalMatrix[Abs[e]^r Sign[e]].v).p),{i,1,dim}]+Table[Total[1/2 y^2  \[Lambda]d Table[v[[i]].(Bdot[[j]].v[[i]]),{i,1,dim}]],{j,1,dim}]//N];


get\[CapitalSigma][q_,r_,ddU_]:=Module[{\[CapitalSigma]=Apply[ddU,q],ve,e,s},
ve=Eigenvectors[\[CapitalSigma]];
e=Eigenvalues[\[CapitalSigma]];
s=Sign[e];
Transpose[ve].DiagonalMatrix[s Abs[e]^r].ve];


HessianH[f_,x_List?VectorQ]:=N[D[f,{x,2}]];
GradientG[f_,x_List?VectorQ]:=N[D[f,{x,1}]];
D3[f_,x_List?VectorQ]:=D[f,{x,3}];


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



