(* ::Package:: *)

BeginPackage["Sampler2`"]

hmc::usage = "hmc[U,Uq,Uqq,Uqqq,Dim,BURNIN,EPISODE,levels,qinit]";

Kp::usage="\!\(\*SubscriptBox[\(K\), \(p\)]\)";
Kq::usage="\!\(\*SubscriptBox[\(K\), \(q\)]\)";
K::usage="K[p,q,i]";

Uq::usage="\!\(\*SubscriptBox[\(U\), \(q\)]\)";
Uqq::usage="\!\(\*SubscriptBox[\(U\), \(qq\)]\)";
Uqqq::usage="\!\(\*SubscriptBox[\(U\), \(qqq\)]\)";


D1::usage="D1[f,x]";
D2::usage="D2[f,x]";
D3::usage="D3[f,x]";
psrf::usage="psrf[QS]";
CHAINS::usage = "mcmc chains, or particles";
dt0::usage = "initial time step";
STEPS::usage="simulation steps"
INTERVAL::usage="INTERVAL for print"
outbnd::usage="check for out of bound"
RATIODT::usage="updating ratio";
RATIOENERGY::usage="updating ratio";
LOWLEVEL::usage="lowest allowable average acceptance probability";
HIGHLEVEL::usage="highest allowable average acceptance probability";

Begin["`Private`"]
CHAINS=3;
STEPS=3;
dt0=0.0000000001;
INTERVAL=1001;
outbnd[q_]:=False;
RATIODT=.1;
RATIOENERGY=.1;
LOWLEVEL=0.1;
HIGHLEVEL=0.9;

hmc[U_,Uq_,Uqq_,Uqqq_,Dim_,BURNIN_,EPISODE_,levels0_,qinit_]:=Module[{
qAll,pAll,UE,ES,Utotal,Ktotal,hi,lo,Htotal,step,Htotals,s,S,AS,u0,anybad,KtotalNew,dt,p,q0,\[Alpha],q,j,bad,i,q1,decayenergy=RATIOENERGY,decaydt=RATIODT,ii=1, QS=ConstantArray[0,{CHAINS( EPISODE-BURNIN),Dim}],ACS={},dts=Table[dt0,Length[levels0]],level,es,h,levels=levels0},
If[qinit!={},
qAll=qinit,
qAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}]];

Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];
Htotals=Table[If[Utotal==0,1,2Utotal],Length[levels]];
level=1;
For[j=1,j\[LessSlantEqual]EPISODE,j++,
pAll=RandomVariate[NormalDistribution[0,1],{CHAINS,Dim}];
KtotalNew = Sum[Apply[K,{pAll[[i]],qAll[[i]],levels[[level]]}],{i,1,CHAINS}];
Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}];

Htotal=Htotals[[level]];
dt=dts[[level]];
Ktotal=Htotal-Utotal;
pAll=pAll Sqrt[Abs[Ktotal/KtotalNew]];

AS={};ES={};
anybad=False;
For[i=1,i<=CHAINS,i++,


bad=False;
p=pAll[[i]];
q=qAll[[i]];
q0=q;

UE={Apply[U,q]};

p=p- dt .5 (Apply[Uq,q]+Apply[Kq,{p,q,levels[[level]]}]);

For[step=1,step<=STEPS && Not[bad],step++,
q1=q;


q=q+dt Apply[Kp,{p,q,levels[[level]]}];

p=p-dt (Apply[Uq,q]+Apply[Kq,{p,q,levels[[level]]}]);

If[outbnd[q],bad=True];

u0=Apply[U,q];
If[Head[u0]!=Real,bad=True];
UE=Append[UE,u0]];

anybad=anybad||bad;

ES=Append[ES,UE];

\[Alpha]=If[bad,0,Exp[Clip[Apply[U,q0]-Apply[U,q],{-10,0}]]];
AS=Append[AS,N[\[Alpha]]];
If[\[Alpha]<RandomVariate[UniformDistribution[]],q=q0];
qAll[[i]]=q;

If[j>BURNIN,
QS[[ii,;;]]=q;
ii=ii+1]];

If[j<=BURNIN && Not[anybad],
s=Union[Flatten[Table[Ordering[ES[[i]],1],{i,1,CHAINS}]]];
S=Union[Flatten[Table[Ordering[ES[[i]],-1],{i,1,CHAINS}]]];
If[s=={1,STEPS+1}&&S=={1,STEPS+1},dt=dt(1+decaydt)];
If[s=={1}&&S=={STEPS+1},dt=dt/(1+decaydt)];
hi=Mean[AS]>HIGHLEVEL;
lo=Mean[AS]<LOWLEVEL;
If[hi||lo,Utotal=Sum[Apply[ U,qAll[[i]]],{i,1,CHAINS}]];
If[hi,Htotal=(Htotal-Utotal)(1+decayenergy)+Utotal,If[lo,Htotal=(Htotal-Utotal)/(1+decayenergy)+Utotal]];
dts[[level]]=dt;
Htotals[[level]]=Htotal;
];
If[Mod[j,INTERVAL]==0,Print[Row[{j,Ktotal,KtotalNew,Mean[AS],StandardDeviation[AS],Utotal,level,s,S,dts,Htotals,dt,lo,hi},"     "]]];
level=Mod[level,Length[levels]]+1];QS];


(*https://mathoverflow.net/questions/229425/derivative-of-eigenvectors-of-a-matrix-with-respect-to-its-components*)
Kq0[p_,q_,r0_,Uqq_,Uqqq_]:=Module[{r=-r0,\[Lambda]d,dim,eig,v,e,y,Bdot,vv,dudq},
\[Lambda]d=r Abs[e]^(r-1);
dim=Length[p];
eig=Eigensystem[Apply[Uqq,q]];
v=Normalize/@eig[[2]];
e=eig[[1]];
y=v.p;
Bdot=Apply[Uqqq,q];vv=Table[Table[Sum[If[i!=j,1/(e[[i]]-e[[j]]) (v[[j]].(Bdot[[k]].v[[i]]))v[[j]],0],{j,1,dim}],{i,1,dim}],{k,1,dim}];dudq=Table[p.(Transpose[vv[[i]]].(DiagonalMatrix[Abs[e]^r Sign[e]].v).p),{i,1,dim}]+Table[Total[1/2 y^2  \[Lambda]d Table[v[[i]].(Bdot[[j]].v[[i]]),{i,1,dim}]],{j,1,dim}]//N];


K[p_,q_,r_]:=Module[{\[CapitalSigma]},
\[CapitalSigma]=getW[q,r,Uqq];
1/2 p.Dot[\[CapitalSigma].p]];

Kp[p_,q_,r_]:=Module[{\[CapitalSigma]},
\[CapitalSigma]=getW[q,r,Uqq];
\[CapitalSigma].p];

Kq[p_,q_,r_]:=Kq0[p,q,r,Uqq,Uqqq];


getW[q_,r_,Uqq_]:=Module[{\[CapitalSigma]=Apply[Uqq,q],ve,e,s,eig},
(*eig=Eigensystem[\[CapitalSigma]];*)
ve=Eigenvectors[\[CapitalSigma]];
e=Eigenvalues[\[CapitalSigma]];
s=Sign[e];
Transpose[ve].DiagonalMatrix[s Abs[e]^(-r)].ve];


D2[f_,x_List?VectorQ]:=N[D[f,{x,2}]];
D1[f_,x_List?VectorQ]:=N[D[f,{x,1}]];
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



