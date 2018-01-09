(* ::Package:: *)

SetDirectory[NotebookDirectory[]];
<<"SingleQubitRB`"


state = vec[DiagonalMatrix[{1,0}]];
meas = Transpose[vec[DiagonalMatrix[{1,0}]]];


infs = 10^RandomReal[{-2,-5},100];
seqLengths = 2^Range[2,11];


(* ::Input::Initialization:: *)
data = Table[
{p,noisyG,delta1,delta2, L} =randomnoisya4[inf];
simData=simRB[noisyG, seqLengths, 500, Inverse[L].state, meas.L];
CI =(1-fitRB[simData, seqLengths]["ParameterConfidenceIntervals"][[1]])/2;
estimatedInfidelity=Mean[CI];
estimatedError=(Max[CI]-Min[CI])/2;
{inf,(1-p[[1]])/2,estimatedInfidelity,estimatedError, delta1,delta2},{inf,infs}];//AbsoluteTiming
Export["data.csv", {{"generatorInfidelity", "trueInfidelity", "estimatedInfidelity","estimatedError", "delta1", "delta2"}}~Join~data, "TextDelimiters"->{"",""}]



