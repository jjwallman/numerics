(* ::Package:: *)

BeginPackage["SingleQubitRB`"]


$DistributedContexts := {$Context, "SingleQubitRB`"}


paulis = {{{1,0},{0,1}},{{0,1},{1,0}},{{0,-I},{I,0}},{{1,0},{0,-1}}};


(* constructing ideal unitaries, t: X->Y->Z *)
T = {{1,-I},{1,I}}/Sqrt[2];
Tmatrices = {{{1,0},{0,1}},T,T.T};
Tperms = {{1,2,3,4},{1,3,4,2},{1,4,2,3}};
pauliproducts = Table[Piecewise[{{a,b==1},{b,a==1},{1,a==b}},24/(a b)],{a,1,4},{b,1,4}];


(* ::Input::Initialization:: *)
randomunitary::usage = "randomunitary[d] returns a Haar-random {d,d}-dimensional unitary matrix based on code from http://www.dr-qubit.org/matlab.html.";


randomsmallunitary::usage = "randomsmallunitary[r] returns a single-qubit unitary with infidelity r to the identity drawn from a unitarily invariant distribution.";


vec::usage = "vec[M] returns the {4,1}-dimensional vector representation of the operator M (typically M is a quantum state) relative to the Pauli basis.";


vecinv::usage = "vecinv[vect] returns the operator with vector representation vect relative to the pauli basis.";


(* ::Input::Initialization:: *)
super::usage = "super[M] returns the {4,4}-dimensional linear representation of the linear map corresponding to conjugation by M (typically M is a unitary or Kraus operator) relative to the Pauli basis.";


rbsequence::usage = "rbsequence[m] returns a single-qubit RB sequence with m random gates.";


randomnoisya4::usage = "randomnoisya4[r] creates a noisy implementation of ideala4 where every Pauli and T matrix is multiplied by an independent sample from randomsmallunitary[r] and then transforms it into the gauge where L = I.";


rgauge::usage = "rgauge[noisygates] puts the noisy implementation of ideala4 into a gauge such that L = I and returns {p, noisygates, delta1, delta2, L} where p is the RB decay rate, delta1 and delta2 characterize the gate-dependence and L is the initial gauge transformation (used to transform states and measurements).";


simRB::usage = "simRB[noisygates, seqlengths, seqsperlength, state, meas] simulates an RB experiment. state and meas should be transformed into the correct gauge";


fitRB::usage = "fitRB[data,seqlengths] fits the output of simRB.";


plotRB::usage = "plotRB[data,seqlengths] plots the output of simRB or simDRB.";


Begin["`Private`"]


randomunitary[d_] := Module[{seeds=QRDecomposition[ArrayReshape[RandomVariate[NormalDistribution[0,1],d^2]+I RandomVariate[NormalDistribution[0,1],d^2],{d,d}]]},
seeds[[2]] = Re[Diagonal[seeds[[2]]]/Abs[Diagonal[seeds[[2]]]]]; 
seeds[[1]].DiagonalMatrix[seeds[[2]]]]


randomsmallunitary[r_] := Module[{V = randomunitary[2],\[Theta]=ArcSin[Sqrt[3 r/2]]},Cos[\[Theta]]{{1,0},{0,1}}-I Sin[\[Theta]]V.{{1,0},{0,-1}}.V\[ConjugateTranspose]]


vec[operator_] := {LinearSolve[Flatten[paulis,{2,3}],Flatten[operator]]}\[Transpose];


vecinv[vect_] := Simplify[Sum[vect[[i,1]]paulis[[i]],{i,4}]];


super[M_] := SparseArray[Sum[vec[M.paulis[[j]].M\[ConjugateTranspose]].SparseArray[{{1,j}->1},{1,4}],{j,1,4}]];


(* ::Input::Initialization:: *)
rbsequence[m_] := Module[{
paulitwirl = RandomInteger[{1,4},m],
Ttwirl = RandomInteger[{0,2},m],
Tcompiled,
paulicompiled},
Tcompiled = Mod[Ttwirl~Join~{0} - {0}~Join~Ttwirl,3];
paulicompiled = Extract[pauliproducts,Transpose[{paulitwirl~Join~{1},{1}~Join~paulitwirl}]];
paulicompiled = Extract[Tperms, Transpose[{1+Mod[-Ttwirl~Join~{0},3],paulicompiled}]];
{1 + Tcompiled,paulicompiled}
]


ideala4 = SparseArray[Table[super[a.b], {a, Tmatrices},{b, paulis}]];
(* ideala4 is a subgroup of the Clifford group that forms a unitary 2-design *)


rgauge[noisygates_] := Module[{p, L, R, ideala42 = ideala4, noisyg, avnoise,deltas},
ideala42[[;;,;;,1,1]] = 0;
{p,L} = Eigensystem[Sum[KroneckerProduct[ideala42[[a,b]],noisygates[[a,b]]],{a,1,3},{b,1,4}]/12,1];
L = Transpose[ArrayReshape[L,{4,4}]] / L[[1,6]];
L[[1,1]] = 1;
noisyg = Table[Inverse[L].noisygates[[a,b]].L, {a,1,3},{b,1,4}];
avnoise = Sum[ConjugateTranspose[ideala4[[a,b]]].noisyg[[a,b]],{a,1,3},{b,1,4}]/12;
deltas = Flatten[Table[Norm[noisyg[[a,b]] - ideala4[[a,b]].avnoise,"Frobenius"],{a,1,3},{b,1,4}]];
{p, noisyg, Max[deltas], Mean[deltas], L}//Chop]


randomnoisya4[r_] := Module[{noisypaulis = Table[b.randomsmallunitary[r],{b,paulis}],
noisyTmatrices = Table[a.randomsmallunitary[r],{a,Tmatrices}]},
rgauge[Table[super[a.b], {a,noisyTmatrices},{b,noisypaulis}]]//Chop]


simRB[noisygates_, seqlengths_, seqsperlength_, state_, meas_] := ParallelTable[2Tr[Fold[Dot, meas, Extract[noisygates, Transpose[ rbsequence[m]]]].state],{m,seqlengths},{k,1,seqsperlength}]


fitRB[data_, seqlengths_]:= Module[{A,B,p,x},
NonlinearModelFit[Transpose[{seqlengths,Mean/@data}],{A + B p^x,0<=A<=1,0<=B<=1,0<=p<=1},{p,{A,0.5},{B,0.5}},x]]


plotRB[data_,seqlengths_,nlm_]:= Module[{means = Mean/@data, stds = (StandardDeviation/@data) / Sqrt[Length/@data]},
ListPlot[(Transpose/@{{seqlengths,means- stds}, {seqlengths,means + stds}, {seqlengths,means}})~Join~{Table[{x,nlm[x]}, {x,seqlengths}]},Filling->{1->{2}},Joined->{False,False,False,True},
PlotStyle->{Gray,Gray,Black,Black},
PlotMarkers->{"-","-","\[FilledCircle]","-"},
Frame->True,
FrameLabel->{"Sequence length", "Survival probability"}]]


End[]


EndPackage[]
