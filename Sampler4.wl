(* ::Package:: *)

BeginPackage["Sampler4`"]

hmc::usage = "hmc3[U,dU,K,dKdp,dKdq,Dim,BURNIN,EPISODE,qinit,levels0]";

dKdp::usage="dKdp[p,q,i]";
dKdq::usage="dKdq[p,q,i]";
K::usage="K[p,q,i]";

U::usage="U[p,q,i]";
dU::usage="dU[p,q,i]";
ddU::usage="ddU[p,q,i]";
dddU::usage="dddU[p,q,i]";

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
REGU::usage="REGU";

Begin["`Private`"]
CHAINS=3;
STEPS=5;
dt0=10^-1;
INTERVAL=1001;
outbnd[q_]:=False;
DECAYDT=.1;
DECAYENERGY=.1;
LOWLEVEL=0.1;
HIGHLEVEL=0.9;
REGU=False;

hmc[Dim_,BURNIN_,EPISODE_,qinit_]:=Module[{
qAll,pAll,UE,ES,Utotal,Ktotal,hi,lo,Htotal,s,S,AS,KtotalNew,dt=dt0,p,q0,\[Alpha],q,j,bad,switch1,k0,h0,k1,h1,i,q1,decayenergy=DECAYENERGY,decaydt=DECAYDT,ii=1, QS=ConstantArray[0,{CHAINS( EPISODE-BURNIN),Dim}],ACS={},es,h},
If[qinit!={},
qAll=qinit,
qAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}]];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotal=If[Utotal==0,1,2Utotal];
For[j=1,j\[LessSlantEqual]EPISODE,j++,
pAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}];
KtotalNew = Sum[Apply[K,{pAll[[i]],qAll[[i]]}],{i,1,CHAINS}];
Utotal=Sum[Apply[U,qAll[[i]]],{i,1,CHAINS}];
Ktotal=Htotal-Utotal;
pAll=pAll Sqrt[Abs[Ktotal/KtotalNew]];
AS={};ES={};
For[i=1,i<=CHAINS,i++,
bad=False;
p=pAll[[i]];
q=qAll[[i]];
q0=q;
If[j<=BURNIN,UE={Apply[U,q]}];

If[REGU,k0=Apply[K,{p,q}];h0=Apply[U,q]+k0,p=p- dt .5 (Apply[dU,q]+Apply[dKdq,{p,q}])];

Do[q1=q;
q=q+dt Apply[dKdp,{p,q}];
p=p-dt (Apply[dU,q]+Apply[dKdq,{p,q}]);

If[REGU,k1=Apply[K,{p,q}];h1=Apply[U,q]+k1;p=p Sqrt[Abs[(k1-(h1-h0))/k1]]];

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
If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt(1+decaydt)];
If[s=={1}&&S=={STEPS+1},dt=dt/(1+decaydt)];
hi=Quantile[AS,.5]>HIGHLEVEL;
lo=Quantile[AS,.5]<LOWLEVEL;
If[hi||lo,Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}]];
If[hi,Htotal=(Htotal-Utotal)(1+decayenergy)+Utotal,If[lo,Htotal=(Htotal-Utotal)/(1+decayenergy)+Utotal]];
];
If[Mod[j,INTERVAL]==0,Print[Row[{j,Ktotal,KtotalNew,Mean[AS],StandardDeviation[AS],Utotal,s,S,dt,Htotal,dt,hi,lo},"     "]]]];QS];



(*https://mathoverflow.net/questions/229425/derivative-of-eigenvectors-of-a-matrix-with-respect-to-its-components*)
dKdq0[p_,q_,ddU_,dddU_]:=Module[{dim,eig,v,e,y,Bdot,vv,dudq,f,fd},
dim=Length[p];
eig=Eigensystem[Apply[ddU,q]];
v=Normalize/@eig[[2]];
e=eig[[1]];
(*f=1/e/Sqrt[Abs[e]];
fd = -3/2 Abs[e]^(-5/2);*)
f=Sign[e]/Sqrt[Abs[e]];
fd=-1/2 Abs[e]^(-3/2);
(*fd=-1/e^2;*)

y=v.p;
Bdot=Apply[dddU,q];
vv=Table[Table[Sum[If[i!=j,1/(e[[i]]-e[[j]]) (v[[j]].(Bdot[[k]].v[[i]]))v[[j]],0],{j,1,dim}],{i,1,dim}],{k,1,dim}];
Table[p.(Transpose[vv[[j]]].(DiagonalMatrix[f].v).p),{j,1,dim}]+Table[Total[1/2 y^2  fd Table[v[[i]].(Bdot[[j]].v[[i]]),{i,1,dim}]],{j,1,dim}]];

K[p_,q_]:=Module[{\[CapitalSigma]},
\[CapitalSigma]=get\[CapitalSigma][q,ddU];
1/2 p.Dot[\[CapitalSigma].p]];

dKdp[p_,q_]:=Module[{\[CapitalSigma]},
\[CapitalSigma]=get\[CapitalSigma][q,ddU];
\[CapitalSigma].p];

dKdq[p_,q_]:=dKdq0[p,q,ddU,dddU];


get\[CapitalSigma][q_,ddU_]:=Module[{\[CapitalSigma]=Apply[ddU,q],ve,e,s},
ve=Eigenvectors[\[CapitalSigma]];
e=Eigenvalues[\[CapitalSigma]];

(*s=1/e;*)
Transpose[ve].DiagonalMatrix[Sign[e]/Sqrt[Abs[e]]].ve];


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



