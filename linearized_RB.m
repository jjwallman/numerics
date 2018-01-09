(* ::Package:: *)

SetDirectory[NotebookDirectory[]];
<<"SingleQubitRB`"


LSestimator[q1_,q2_,m1_,m2_] := {(q1^m2/q2^m1)^(1/(m2-m1)), (q2/q1)^(1/(m2-m1))}


(* ::Input::Initialization:: *)
binp[x_, n_, p_]:=Binomial[n,x]p^x (1-p)^(n-x)


(* ::Input::Initialization:: *)
estimators[m_, k_] := estimators[m, k] =Transpose[Sort[Flatten[Table[{N[((2y -k)/(2x-k))^(1/(m-4))],x,y},{x,Range[Ceiling[(k+1)/2],k]},{y,Range[Ceiling[(k+1)/2],k]}],1],#1[[1]]<#2[[1]]&]]


(* ::Input::Initialization:: *)
cdf[p_,A_,k_] := Module[{m = Ceiling[1/(2-2p)],ests,x,y,probs},
{ests,x,y}= estimators[m,k];
probs =  Thread[binp[x,k,A p^4 + .5]]Thread[binp[y,k,A p^m + .5]];
Transpose[{ests,FoldList[Plus,probs]}]]


(* ::Input::Initialization:: *)
exportcdf[p_,A_,k_,title_] :=
Module[{data = cdf[p,A,k],cil,ciu},
data[[;;,1]]=(1-data[[;;,1]])/(1-p);
cil=SelectFirst[data,#[[2]]>=0.025&][[1]];
ciu=SelectFirst[data,#[[2]]>=0.975&][[1]];
Export["cdfp"<>StringReplace[ToString[p],"."->"_"]<>"A"<>StringReplace[ToString[A],"."->"_"]<>"samples"<>ToString[k]<>".pdf",ListPlot[data,
PlotRange->{{0,4},{0,1}},
FrameLabel->{{ToExpression["\\mathrm{Pr}(\\hat{p} \\geq x)", TeXForm,HoldForm],None},{ToExpression["\\frac{1-x}{1-p}", TeXForm,HoldForm],Row[{title[[1]], Spacer@title[[2]]}]}},
Ticks->{Automatic,ScientificForm},
AxesStyle->Black,
Frame->True,
RotateLabel->False,
Epilog->{{Gray,Line@{{cil,0},{cil,1}}},{Gray,Line@{{1,0},{1,1}}},{Gray,Line@{{ciu,0},{ciu,1}}}}]]]


(* ::Code:: *)
(*labels = <|{0.35, 0.99}->{"a) A=0.35, p=0.99", 280},*)
(*{0.49, 0.99}->{"b) A=0.49, p=0.99", 280},*)
(*{0.35, 0.9999}->{"c) A=0.35, p=0.9999", 270},*)
(*{0.49, 0.9999}->{"d) A=0.49, p=0.9999", 270}|>;*)
(*Do[exportcdf[p,A,100,labels[{A,p}]], {p,{0.99,0.9999}},{A,{0.35,0.49}}];*)


(* ::Code:: *)
(*CIwidth = Table[Table[temp=cdf[p,.49,numseqs];*)
(*cil=SelectFirst[temp,#[[2]]>=0.025&][[1]];*)
(*ciu=SelectFirst[temp,#[[2]]>=0.975&][[1]];*)
(*{numseqs,(ciu-cil)/(1-p)},*)
(*{numseqs,Round[10^RandomReal[{1.3,3},20]]}],*)
(*{p,{0.99,.99999999}}];*)
(*Export["CIplot.pdf",ListLogLinearPlot[CIwidth,*)
(*AxesStyle->Black,*)
(*PlotMarkers->{"+","\[FilledCircle]"},*)
(*Frame->True,*)
(*RotateLabel->False,*)
(*FrameLabel->{"k",ToExpression["\\frac{\\delta}{1-p}", TeXForm,HoldForm]}]]*)


posterior[q1_,q2_,m1_,m2_,k_] := Module[{hatA,hatp,post},
{hatA,hatp} = LSestimator[q1,q2,m1,m2];
post=Table[{p,NIntegrate[(.5+A p^m1)^(q1 k) (.5- A p^m1)^(k - q1 k) (.5+A p^m2)^(q2 k) (.5- A p^m2)^(k - q2 k),{A,0,0.5}]},{p,Subdivide[10 hatp - 9,1,100]}];
post = Transpose[{1-post[[;;,1]],FoldList[Plus,post[[;;,2]]]/Total[post[[;;,2]]]}]]


Export["posterior.pdf",ListPlot[{posterior[.98,.7,4,1000,50],posterior[.98,.7,4,1000,500]},
AxesStyle->Black,
Joined->True,
PlotStyle -> {Automatic, Dashed},
Frame->True,
RotateLabel->False,
FrameLabel->{"1-x",ToExpression["\\mathrm{Pr}(\\hat{p}\\geq x|q_1,q_2)", TeXForm,HoldForm]},
Epilog->{Inset["500 samples", {0.0005, 0.1}],Inset["50 samples", {0.003, 0.1}]}]]


(* ::Code:: *)
(*state = vec[DiagonalMatrix[{1,0}]];*)
(*meas = Transpose[vec[DiagonalMatrix[{1,0}]]];*)
(*inf = 1. 10^-3;*)
(*nshots = 100;*)
(*nsamples = 30;*)
(*CIs = Table[*)
(*{p,noisyG,delta1,delta2, L} =randomnoisya4[inf];*)
(*p = p[[1]];*)
(*m = Round[1/(2-2p)];*)
(*Lstate = Inverse[L].state;*)
(*Lmeas = meas.L;*)
(*data = simDRB[noisyG, {4, m}, nsamples, Lstate, Lmeas];*)
(*empiricalweights = Map[RandomVariate[BinomialDistribution[nshots,#]]&,data,{2}]/nshots;*)
(*bootstrapdist = Table[EmpiricalDistribution[j],{j,empiricalweights}];*)
(*estimates = Mean/@empiricalweights;*)
(*estimates = LSestimator[estimates[[1]],estimates[[2]],4,m];*)
(*bootstrappars = Table[LSestimator[N[Mean[RandomVariate[bootstrapdist[[1]],nsamples]]]-0.5,N[Mean[RandomVariate[bootstrapdist[[2]],nsamples]]]-0.5,4,m],{j,1,2000}];*)
(*bootp = Sort[bootstrappars[[;;,2]]];*)
(*{p}~Join~bootp[[{Floor[.025 * Length[bootp]],Ceiling[.975 * Length[bootp]]}]],{run,1,50}];*)
(*CIs = SortBy[CIs,First];*)


(* ::Code:: *)
(*CIplot=ListPlot[{CIs[[;;,1;;2]],CIs[[;;,{1,3}]]},*)
(*Joined->{False,False},*)
(*Filling->{1->{2}},*)
(*PlotStyle->Black,*)
(*PlotMarkers -> {"-", "-"},*)
(*Epilog->{Line[{{0,0},{1,1}}]},*)
(*AxesStyle->Black,*)
(*Frame->True]*)
(*Export["bootstrap_plot.pdf",CIplot]*)



