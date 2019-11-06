(* ::Package:: *)

(* ::Section:: *)
(*HazmaTools*)


Global`$LoadFeynArts=True;
Needs["FeynCalc`"]
$FAVerbose=False;


BeginPackage["HazmaTools`"];


$HazmaModel::usage="$HazmaModel denotes the model which should be used. Set `$HazmaModel=\"scalar\"` for scalar mediator, `$HazmaModel=\"vector\"` for vector mediator, etc. By default, it is set to scalar."


state\[Pi]0::usage="FeynArts state for neutral pion";
state\[Pi]m::usage="FeynArts state for negatively charged-pion";
state\[Pi]p::usage="FeynArts state for positively charged-pion";
state\[Nu]::usage="FeynArts state for neutrino"; 
state\[Nu]bar::usage="FeynArts state for anti-neutrino"; 
statel::usage="FeynArts state for down-type lepton"; 
statelbar::usage="FeynArts state for anit down-type lepton";
statek0::usage="FeynArts state for neutral kaon:";
statek0bar::usage="FeynArts state for anti neutral kaon"; 
statekm::usage="FeynArts state for negatively charged kaon"; 
statekp::usage="FeynArts state for positively charged kaon";
state\[Eta]::usage="FeynArts state for eta meson";
state\[Rho]::usage="FeynArts state for rho vector meson";
statex::usage="FeynArts state for dark matter";
statexbar::usage="FeynArts state for anit dark matter";
states::usage="FeynArts state for scalar mediator";
state\[Gamma]::usage="FeynArts state for photon"
statep::usage="FeynArts state for psuedo-scalar mediator";
statev::usage="FeynArts state for vector mediator";


HazmaGenerateDiagrams::usage="HazmaGenerateDiagrams[inStates,outStates] Generate `FeynArts` Feynmann diagrams for the process `inStates`->`outStates`.";


HazmaSubstitutionUnstableMediatorPropagator::usage="HazmaSubstitutionUnstableMediatorPropagator[exp] Substitue the mediator propagator for its unstable form in `exp`";


HazmaCreateAmplitude::usage="HazmaCreateAmplitude[inStates,outStates] creates an amplitude for `inStates`->`outStates` in the Scalar-Mediator model.";


HazmaCreateAmplitudeSquared::usage="HazmaCreateAmplitudeSquared[inStates,outStates] creates the squared spin-averaged amplitude for `inStates`->`outStates` in the Scalar-Mediator model.";


HazmaComputeCrossSection22::usage="HazmaComputeCrossSection22[inStates,outStates,Q] computes 2->2 cross section forprocess `inStates`->`outStates` given a center of mass energy `Q` in the Scalar-Mediator model.";


HazmaComputeDNDE::usage="HazmaComputeDNDE[inStates,outStates,Q] computes \[Gamma]-ray spectrum for `inStates`->`outStates` + \[Gamma] given a center of mass energy `Q` in the Scalar-Mediator model. Returns an expression in terms of ";


SubstituteMasses::usage"If true, masses are subsituted, otherwise the are left as scalar products of momenta."


Begin["`Private`"]


Print["MeV-DM - Scalar Mediator: version 0.1"];
Print["by Adam Coogan and Logan A. Morrison"];


$HazmaModel="scalar";


(* ::Subsection::Closed:: *)
(*States*)


(* SM *)
state\[Pi]0=FeynArts`S[1];
state\[Pi]m=FeynArts`S[2];
state\[Pi]p=-FeynArts`S[2];
statek0=FeynArts`S[3]; 
statek0bar=-FeynArts`S[3]; 
statekm=FeynArts`S[4]; 
statekp=-FeynArts`S[4];
state\[Eta]=FeynArts`S[5];
state\[Gamma]=FeynArts`V[1];
state\[Rho]=FeynArts`V[2];
state\[Omega]=FeynArts`V[3];
(* Leptons *)
state\[Nu]=FeynArts`F[1]; 
state\[Nu]bar=-FeynArts`F[1]; 
statel=FeynArts`F[2]; 
statelbar=-FeynArts`F[2];
(* DM *)
statex=FeynArts`F[3]; 
statexbar=-FeynArts`F[3];
(* mediators *)
states=FeynArts`S[7];
statep=FeynArts`S[7];
statev=FeynArts`V[4];


(* SM *)
FeynArts`Mass[state\[Pi]0]=FeynArts`TheMass[state\[Pi]0];
FeynArts`Mass[state\[Pi]m]=FeynArts`TheMass[state\[Pi]m];
FeynArts`Mass[state\[Pi]p]=FeynArts`TheMass[state\[Pi]p];
FeynArts`Mass[statek0]=FeynArts`TheMass[statek0]; 
FeynArts`Mass[statek0bar]=FeynArts`TheMass[statek0bar]; 
FeynArts`Mass[statekm]=FeynArts`TheMass[statekm]; 
FeynArts`Mass[statekp]=FeynArts`TheMass[statekp];
FeynArts`Mass[state\[Eta]]=FeynArts`TheMass[state\[Eta]];
FeynArts`Mass[state\[Gamma]]=FeynArts`TheMass[state\[Gamma]];
FeynArts`Mass[state\[Rho]]=FeynArts`TheMass[state\[Rho]];
FeynArts`Mass[state\[Omega]]=FeynArts`TheMass[state\[Omega]];
(* Leptons *)
FeynArts`Mass[state\[Nu]]=FeynArts`TheMass[state\[Nu]]; 
FeynArts`Mass[state\[Nu]bar]=FeynArts`TheMass[state\[Nu]bar]; 
FeynArts`Mass[statel]=FeynArts`TheMass[statel]; 
FeynArts`Mass[statelbar]=FeynArts`TheMass[statelbar];
(* DM *)
FeynArts`Mass[statex]=FeynArts`TheMass[statex]; 
FeynArts`Mass[statexbar]=FeynArts`TheMass[statexbar];
(* mediators *)
FeynArts`Mass[states]=FeynArts`TheMass[states];
FeynArts`Mass[statep]=FeynArts`TheMass[statep];
FeynArts`Mass[statev]=FeynArts`TheMass[statev];


(* determine if state is a fermion *)
isFermion[state_]:=Module[{head},
	(* if state looks like -F[1], need to strip off the `Times` head *)
	head=If[Head[state]===Times, Head[ReplaceAll[state,Times[a_,b__]:>b]], Head[state]];
	head===FeynArts`F
]
(* determine if state is a vector *)
isVector[state_]:=Module[{head},
	(* if state looks like -V[1], need to strip off the `Times` head *)
	head=If[Head[state]===Times, Head[ReplaceAll[state,Times[a_,b__]:>b]], Head[state]];
	head===FeynArts`V
]	


(* ::Subsection::Closed:: *)
(*Generating Feynman Diagrams from FeynArts*)


GetModelName[]:=Module[{modelName},
	If[$HazmaModel==="scalar",Return[{"EFT_MeV_DM_scalar"<>"/"<>"EFT_MeV_DM_scalar"}]];
	If[$HazmaModel==="vector",Return[{"EFT_MeV_DM_vector_no_contact"<>"/"<>"EFT_MeV_DM_vector_no_contact"}]];
	Throw[$Failed, Abort::InvalidHazmaModel]
]


GetGenericModelName[]:=Module[{modelName},
	If[$HazmaModel==="scalar",Return[{"Lorentz","EFT_MeV_DM_scalar"<>"/"<>"EFT_MeV_DM_scalar"}]];
	If[$HazmaModel==="vector",Return[{"Lorentz","EFT_MeV_DM_vector_no_contact"<>"/"<>"EFT_MeV_DM_vector_no_contact"}]];
	Throw[$Failed, Abort::InvalidHazmaModel]
]


Options[HazmaGenerateDiagrams]={
	LoopNumber->0,
	Adjacencies->{3,4,5},
	ExcludeParticles->{},
	ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	InsertionLevel->{FeynArts`Particles},
	Paint->False
};

HazmaGenerateDiagrams[inStates_, outStates_,OptionsPattern[]]:=Block[{tops,diags},
	tops = FeynArts`CreateTopologies[
		OptionValue[LoopNumber], 
		Length[inStates]-> Length[outStates],
		FeynArts`ExcludeTopologies-> OptionValue[ExcludeTopologies],
		FeynArts`Adjacencies-> OptionValue[Adjacencies]
	];
	diags=FeynArts`InsertFields[
		tops,
		inStates -> outStates,
		FeynArts`InsertionLevel -> OptionValue[InsertionLevel],
		FeynArts`Model -> GetModelName[],
		FeynArts`GenericModel -> GetGenericModelName[],
		FeynArts`ExcludeParticles -> OptionValue[ExcludeParticles]
	];
	If[OptionValue[Paint],FeynArts`Paint[diags],FeynArts`SheetHeader->None];
	diags
]


(* ::Subsection::Closed:: *)
(*Computing Amplitudes*)


SubstitutionUnstableScalarPropagator[exp_]:=ReplaceAll[exp,{
FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum+b__FeynCalc`Momentum,Global`ms]]:>1/(FeynCalc`SP[a+b,a+b]-Global`ms^2+I Global`widths Global`ms),
FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a_FeynCalc`Momentum-b_FeynCalc`Momentum,Global`ms]]:>1/(FeynCalc`SP[a-b,a-b]-Global`ms^2+I Global`widths Global`ms),
	FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum,Global`ms]]:>1/(FeynCalc`SP[a,a]-Global`ms^2+I Global`widths Global`ms)}];
	
SubstitutionUnstableVectorPropagator[exp_]:=ReplaceAll[exp,{
FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum+b__FeynCalc`Momentum,Global`mv]]:>1/(FeynCalc`SP[a+b,a+b]-Global`mv^2+I Global`widthv Global`mv),
FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a_FeynCalc`Momentum-b_FeynCalc`Momentum,Global`mv]]:>1/(FeynCalc`SP[a-b,a-b]-Global`mv^2+I Global`widthv Global`mv),
	FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum,Global`mv]]:>1/(FeynCalc`SP[a,a]-Global`mv^2+I Global`widthv Global`mv)}];
	
HazmaSubstitutionUnstableMediatorPropagator[exp_]:=Module[{},
	If[$HazmaModel==="scalar",Return[SubstitutionUnstableScalarPropagator[exp]]];
	If[$HazmaModel==="vector",Return[SubstitutionUnstableVectorPropagator[exp]]];
	Throw[$Failed, Abort::InvalidHazmaModel]
]


Clear[HazmaCreateAmplitude];

Options[HazmaCreateAmplitude]={
	LoopNumber->0,
	Adjacencies->{3,4,5},
	ExcludeParticles->{},
	ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	InsertionLevel->{FeynArts`Particles},
	Paint->False,
	ChangeDimension->4,
	FinalSubstitutions->{Global`M$FACouplings},
	IncomingMomenta->{},
	List->False,
	LoopMomenta->{},
	LorentzIndexNames->{Global`\[Mu],Global`\[Nu],Global`\[Alpha],Global`\[Beta],Global`\[Rho],Global`\[Sigma]},
	OutgoingMomenta->{},
	SMP->True,
	TransversePolarizationVectors->{},
	UndoChiralSplittings->True,
	GaugeRules-> _FeynArts`FAGaugeXi->1,
	PreFactor->-I (2 \[Pi])^(-4 * FeynArts`LoopNumber),
	Truncated->False
};

HazmaCreateAmplitude[inStates_,outStates_,OptionsPattern[]] := Module[{tops, ampFA,ampFC, diagrams},

FeynCalc`ClearScalarProducts[];
(* create the topologies based on input states *)
diagrams=HazmaGenerateDiagrams[inStates,outStates,
	LoopNumber->OptionValue[LoopNumber],
	Adjacencies->OptionValue[Adjacencies],
	ExcludeParticles->OptionValue[ExcludeParticles],
	ExcludeTopologies->OptionValue[ExcludeTopologies],
	InsertionLevel->OptionValue[InsertionLevel],
	Paint->OptionValue[Paint]
];
	
(*Compute the amplitude *)
ampFA=FeynArts`CreateFeynAmp[diagrams, 
	FeynArts`GaugeRules->OptionValue[GaugeRules],
	FeynArts`PreFactor->OptionValue[PreFactor],
	FeynArts`Truncated->OptionValue[Truncated]
];

ampFC=FeynCalc`FCFAConvert[
	ampFA,
	FeynCalc`IncomingMomenta-> OptionValue[IncomingMomenta],
	FeynCalc`OutgoingMomenta-> OptionValue[OutgoingMomenta],
	FeynCalc`LorentzIndexNames-> OptionValue[LorentzIndexNames],
	FeynCalc`LoopMomenta->OptionValue[LoopMomenta], 
	FeynCalc`UndoChiralSplittings->OptionValue[UndoChiralSplittings],
	FeynCalc`SMP-> OptionValue[SMP],
	FeynCalc`ChangeDimension-> OptionValue[ChangeDimension],
	FeynCalc`FinalSubstitutions-> OptionValue[FinalSubstitutions],
	FeynCalc`TransversePolarizationVectors->OptionValue[TransversePolarizationVectors]
];

ampFC=Flatten[ampFC];
ampFC=If[OptionValue[List],ampFC,Total[ampFC]];
ampFC=ReplaceAll[ampFC,{FeynArts`FALeviCivita->FeynCalc`Eps,Global`FAFourVector->FeynCalc`Momentum}];
(* subsitute unstable propagator *)
ampFC=Map[HazmaSubstitutionUnstableMediatorPropagator,ampFC];
ampFC=Map[FeynCalc`PropagatorDenominatorExplicit,ampFC];
ampFC=Map[FeynCalc`Contract,ampFC];
ampFC=Map[FeynCalc`DiracSimplify[#,FeynCalc`DiracSubstitute67->True]&,ampFC];
ampFC=Map[FeynCalc`ExpandScalarProduct,ampFC];
Map[Simplify,ampFC]
];


(* ::Subsection:: *)
(*Computing Squared Amplitudes*)


Options[HazmaCreateAmplitudeSquared]={
	LoopNumber->0,
	Adjacencies->{3,4,5},
	ExcludeParticles->{},
	ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	InsertionLevel->{FeynArts`Particles},
	Paint->False,
	ChangeDimension->4,
	FinalSubstitutions->{Global`M$FACouplings},
	IncomingMomenta->{},
	List->False,
	LoopMomenta->{},
	LorentzIndexNames->{Global`\[Mu],Global`\[Nu],Global`\[Alpha],Global`\[Beta],Global`\[Rho],Global`\[Sigma]},
	OutgoingMomenta->{},
	SMP->True,
	TransversePolarizationVectors->{},
	UndoChiralSplittings->True,
	GaugeRules->_FeynArts`FAGaugeXi->1,
	PreFactor->-I (2 \[Pi])^(-4 FeynArts`LoopNumber),
	Truncated->False
};

HazmaCreateAmplitudeSquared[inStates_, outStates_,OptionsPattern[]] := Module[{msqrd, amp,dim,momentum,inMomenta,outMomenta},

(* Create the amplitude *)
amp=HazmaCreateAmplitude[
	inStates, 
	outStates,
	LoopNumber->OptionValue[LoopNumber],
	Adjacencies->OptionValue[Adjacencies],
	ExcludeParticles->OptionValue[ExcludeParticles],
	ExcludeTopologies->OptionValue[ExcludeTopologies],
	InsertionLevel->OptionValue[InsertionLevel],
	Paint->OptionValue[Paint],
	ChangeDimension->OptionValue[ChangeDimension],
	FinalSubstitutions->OptionValue[FinalSubstitutions],
	IncomingMomenta->OptionValue[IncomingMomenta],
	List->OptionValue[List],
	LoopMomenta->OptionValue[LoopMomenta],
	LorentzIndexNames->OptionValue[LorentzIndexNames],
	OutgoingMomenta->OptionValue[OutgoingMomenta],
	SMP->OptionValue[SMP],
	TransversePolarizationVectors->OptionValue[TransversePolarizationVectors],
	UndoChiralSplittings->OptionValue[UndoChiralSplittings],
	GaugeRules->OptionValue[GaugeRules],
	PreFactor->OptionValue[PreFactor],
	Truncated->OptionValue[Truncated]
	];

FeynCalc`ClearScalarProducts[];
(* determine incoming/outgoing momenta *)
inMomenta=OptionValue[IncomingMomenta];
outMomenta=OptionValue[OutgoingMomenta];

(* Fix momenta - If user didn't specify momenta, populate them based on what FeynCalc does. *)
For[i=Length[inMomenta]+1,i<=Length[inStates],i++,
	AppendTo[inMomenta,ToExpression["InMom"<>ToString[i]]]
];
For[i=Length[outMomenta]+1,i<=Length[outStates],i++,
	AppendTo[outMomenta,ToExpression["OutMom"<>ToString[i]]]
];
(* set masses *)
For[i=1,i<=Length[inMomenta],i++,
	FeynCalc`SP[inMomenta[[i]]]=FeynArts`TheMass[inStates[[i]]]^2;
];
For[i=1,i<=Length[outMomenta],i++,
	FeynCalc`SP[outMomenta[[i]]]=FeynArts`TheMass[outStates[[i]]]^2
];

(* Square it *)
msqrd=amp*FeynCalc`FCRenameDummyIndices[FeynCalc`ComplexConjugate[amp]];

(* If any of the in states are fermionic, average over spins *)
For[i=1,i<=Length[inStates],i++,
	If[isFermion[inStates[[i]]], 
		msqrd=FeynCalc`FermionSpinSum[msqrd,FeynCalc`ExtraFactor->1/2],
		Continue
	];
];

(* If any of the out states are fermionic, sum over spins *)
For[i=1,i<=Length[outStates],i++,
	If[isFermion[outStates[[i]]], 
		msqrd=FeynCalc`FermionSpinSum[msqrd],
		Continue
	];
];

(* actually perform the fermion spin sums *)
msqrd=ReplaceAll[msqrd,FeynCalc`DiracTrace->FeynCalc`TR];

(* If any of the states are vectors, do a sum over polarizations. *)
(* Add 1/2 (1/3) massless (massive) inital state vectrors *)
For[i=1,i<=Length[inStates],i++,
	If[isVector[inStates[[i]]], 
		If[FeynArts`TheMass[inStates[[i]]]===0,
			msqrd=(1/2)FeynCalc`DoPolarizationSums[msqrd,inMomenta[[i]],0],
			msqrd=(1/3)FeynCalc`DoPolarizationSums[msqrd,inMomenta[[i]]];
		], 
		Continue
	];
];

For[i=1,i<=Length[outStates],i++,
	If[isVector[outStates[[i]]], 
		If[FeynArts`TheMass[outStates[[i]]]===0,
			msqrd=FeynCalc`DoPolarizationSums[msqrd,outMomenta[[i]],0],
			msqrd=FeynCalc`DoPolarizationSums[msqrd,outMomenta[[i]]];
		], 
		Continue
	];
];

msqrd=FeynCalc`Contract[msqrd];
msqrd=Simplify[msqrd];
FeynCalc`ClearScalarProducts[];
msqrd

];


(* ::Subsection:: *)
(*Compute 2->2 Cross Sections*)


ComputeCrossSectionPreFactor[Q_,m1_,m2_]:=Module[{E1},
	E1=(Q^2+m1^2-m2^2)/(2*Q);
	1/(4*Q*Sqrt[E1^2-m1^2])
]


Options[HazmaComputeCrossSection22]={
	LoopNumber->0,
	Adjacencies->{3,4,5},
	ExcludeParticles->{},
	ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	InsertionLevel->{FeynArts`Particles},
	Paint->False,
	ChangeDimension->4,
	FinalSubstitutions->{Global`M$FACouplings},
	IncomingMomenta->{},
	List->False,
	LoopMomenta->{},
	LorentzIndexNames->{Global`\[Mu],Global`\[Nu],Global`\[Alpha],Global`\[Beta],Global`\[Rho],Global`\[Sigma]},
	OutgoingMomenta->{},
	SMP->True,
	TransversePolarizationVectors->{},
	UndoChiralSplittings->True,
	GaugeRules->_FeynArts`FAGaugeXi->1,
	PreFactor->-I (2 \[Pi])^(-4 FeynArts`LoopNumber),
	Truncated->False
};

HazmaComputeCrossSection22[inStates_, outStates_,Q_,OptionsPattern[]] := Module[{msqrd,inMomenta,outMomenta,m1,m2,m3,m4,s,t,u,E1,E3,z,p1,p3},
	ClearScalarProducts[];
	(* Create the amplitude *)
	msqrd=HazmaCreateAmplitudeSquared[
		inStates, 
		outStates,
		LoopNumber->OptionValue[LoopNumber],
		Adjacencies->OptionValue[Adjacencies],
		ExcludeParticles->OptionValue[ExcludeParticles],
		ExcludeTopologies->OptionValue[ExcludeTopologies],
		InsertionLevel->OptionValue[InsertionLevel],
		Paint->OptionValue[Paint],
		ChangeDimension->OptionValue[ChangeDimension],
		FinalSubstitutions->OptionValue[FinalSubstitutions],
		IncomingMomenta->OptionValue[IncomingMomenta],
		List->OptionValue[List],
		LoopMomenta->OptionValue[LoopMomenta],
		LorentzIndexNames->OptionValue[LorentzIndexNames],
		OutgoingMomenta->OptionValue[OutgoingMomenta],
		SMP->OptionValue[SMP],
		TransversePolarizationVectors->OptionValue[TransversePolarizationVectors],
		UndoChiralSplittings->OptionValue[UndoChiralSplittings],
		GaugeRules->OptionValue[GaugeRules],
		PreFactor->OptionValue[PreFactor],
		Truncated->OptionValue[Truncated]
	];
	(* determine incoming/outgoing momenta *)
	inMomenta=OptionValue[IncomingMomenta];
	outMomenta=OptionValue[OutgoingMomenta];
	(* Fix momenta *)
	For[i=Length[inMomenta]+1,i<=Length[inStates],i++,
		AppendTo[inMomenta,ToExpression["InMom"<>ToString[i]]]
	];
	For[i=Length[outMomenta]+1,i<=Length[outStates],i++,
		AppendTo[outMomenta,ToExpression["OutMom"<>ToString[i]]]
	];
	
	m1=FeynArts`TheMass[inStates[[1]]];
	m2=FeynArts`TheMass[inStates[[2]]];
	m3=FeynArts`TheMass[outStates[[1]]];
	m4=FeynArts`TheMass[outStates[[2]]];

	FeynCalc`SetMandelstam[s,t,u,inMomenta[[1]],inMomenta[[2]],-outMomenta[[1]],-outMomenta[[2]],m1,m2,m3,m4];
	msqrd=Simplify[msqrd];

	msqrd=ReplaceAll[msqrd,{u->m1^2+m2^2+m3^2+m4^2-s-t}];

	E1=(s+m1^2-m2^2)/(2Sqrt[s]);
	E3=(s+m3^2-m4^2)/(2Sqrt[s]);
	p1=Sqrt[E1^2-m1^2];
	p3=Sqrt[E3^2-m3^2];
	msqrd=Simplify[ReplaceAll[msqrd,{t->m1^2+m3^2-2(E1*E3-z p1*p3)}]];

	ReplaceAll[(2\[Pi])/(64\[Pi]^2 Q^2) (p3/p1) Integrate[msqrd,{z,-1,1},GenerateConditions->False],{s->Q^2}]

];


(* ::Subsection:: *)
(*Computing dNdE*)


ComputeMandelstamTBound[M_, m1_, m2_, m3_, s_, t_]:=Module[{E1,E2},
	E1=(M^2+m1^2-s)/(2*M);
	E2=(M^2+m2^2-t)/(2*M);
	Solve[M^2+2*E1*E2+m1^2+m2^2-m3^2-2*M*(E1+E2)==2Sqrt[E1^2-m1^2] Sqrt[E2^2-m2^2],t]
];


Options[ScalarMediatorComputeDNDE]={
	Adjacencies->{3,4,5},
	Paint->False,
	FinalSubstitutions->{Global`M$FACouplings}
};

ScalarMediatorComputeDNDE[inStates_,outStates_,Q_,OptionsPattern[]]:=Module[{newOutStates,msqrd,\[Sigma]0,p1,p2,p3,p4,k,m1,m2,m3,m4,s,t,u,E1,E3,tbounds,tintegral,dndE,preFactor},
	\[Sigma]0=HazmaComputeCrossSection22[inStates,outStates,Q,Adjacencies->OptionValue[Adjacencies]];
	newOutStates=outStates;
	AppendTo[newOutStates,state\[Gamma]];
	m1=FeynArts`TheMass[inStates[[1]]];
	m2=FeynArts`TheMass[inStates[[2]]];
	m3=FeynArts`TheMass[outStates[[1]]];
	m4=FeynArts`TheMass[outStates[[2]]];
	(* set masses *)
	FeynCalc`SP[p1,p1]=m1^2;
	FeynCalc`SP[p2,p2]=m2^2;
	FeynCalc`SP[p3,p3]=m3^2;
	FeynCalc`SP[p4,p4]=m4^2;
	FeynCalc`SP[k,k]=0;
	msqrd=HazmaCreateAmplitudeSquared[
		inStates, 
		newOutStates,
		Adjacencies->OptionValue[Adjacencies],
		Paint->OptionValue[Paint],
		ChangeDimension->4,
		FinalSubstitutions->Global`M$FACouplings,
		IncomingMomenta->{p1,p2},
		List->False,
		LorentzIndexNames->{\[Mu],\[Nu],\[Alpha],\[Beta]},
		OutgoingMomenta->{p3,p4,k},
		SMP->True,
		TransversePolarizationVectors->{k},
		(* ignore diagrams with internal photon. They are higher order in \[Alpha] *)
		ExcludeParticles->{state\[Gamma]}
	];
	(* Set the Mandelstam variables *)
	(* s = (P - k)^2 = (p3 + p4)^2 = m3^2 + m4^2 + 2p3.p4 *)
	FeynCalc`SP[p3,p4]=(s-m3^2-m4^2)/2;
	(* t = (P - p3)^2 = (p4 + k)^2 = m4^2 + 2k.p4 *)
	FeynCalc`SP[k, p4]=(t-m4^2)/2;
	(* u = (P - p4)^2 = (p3 + k)^2 = m3^2 + 2p3.k *)
	FeynCalc`SP[k, p3]=(u-m3^2)/2;

	FeynCalc`SP[p1,p2]=(Q^2-m1^2-m2^2)/2;

	msqrd=msqrd/.{u->Q^2+m3^2+m4^2-s-t};

	tbounds=ComputeMandelstamTBound[Q, 0, m3, m4, s, t];
	E1=(Q^2+m1^2-m2^2)/(2Q);
	tintegral=Integrate[msqrd,t,GenerateConditions->False];
	dndE=Simplify[((tintegral/.tbounds[[1]])-(tintegral/.tbounds[[2]])),Assumptions->{Q>0,Q>m1+m2,Q>m3+m4}];
	preFactor=(2*Q)/(4*Q Sqrt[E1^2-m1^2]) 1/(16*Q^2 (2\[Pi])^3) * 1/\[Sigma]0;
	Simplify[
		ReplaceAll[
			(-1)preFactor*dndE,
			{s->Q^2(1 - 2 Global`E\[Gamma] / Q),Global`qe->Sqrt[4\[Pi] Global`\[Alpha]em]}],
		Assumptions->{Q>0, Q-(m3+m4)^2/Q > 2 * Global`E\[Gamma] > 0,m1>0,m2>0,m3>0,m4>0}]
];


(* ::Text:: *)
(*Vector Mediator*)


Options[VectorMediatorComputeDNDE]={
	Adjacencies->{3,4,5},
	Paint->False,
	FinalSubstitutions->{Global`M$FACouplings}
};

VectorMediatorComputeDNDE[inStates_,outStates_,Q_,OptionsPattern[]]:=Module[{newOutStates,msqrd,\[Sigma]0,P,p1,p2,p3,p4,k,m1,m2,m3,m4,s,t,u,E1,E3,tbounds,tintegral,dndE,preFactor,X\[Mu]\[Nu],L\[Mu]\[Nu],ampFS,ampIS,X,L},
	\[Sigma]0=HazmaComputeCrossSection22[inStates,outStates,Q,Adjacencies->OptionValue[Adjacencies]];
	newOutStates=outStates;
	AppendTo[newOutStates,state\[Gamma]];
	(* Factor the amplitude into L\[Mu]\[Nu]*X\[Mu]\[Nu]/propagator *)
	(* Initial state M(xx\[Rule]V)^2/4 *)
	ampIS=FeynCalc`Contract[HazmaCreateAmplitude[
		inStates, 
		{statev},
		Adjacencies->OptionValue[Adjacencies],
		Paint->OptionValue[Paint],
		ChangeDimension->4,
		FinalSubstitutions->Global`M$FACouplings,
		IncomingMomenta->{p1,p2},
		List->False,
		LorentzIndexNames->{},
		OutgoingMomenta->{P},
		SMP->True,
		TransversePolarizationVectors->{},
		ExcludeParticles->{state\[Gamma]}
	]];
	(* Final state piece *)
	ampFS=FeynCalc`Contract[HazmaCreateAmplitude[
		{statev}, 
		newOutStates,
		Adjacencies->OptionValue[Adjacencies],
		Paint->OptionValue[Paint],
		ChangeDimension->4,
		FinalSubstitutions->Global`M$FACouplings,
		IncomingMomenta->{P},
		List->False,
		LorentzIndexNames->{},
		OutgoingMomenta->{p3,p4,k},
		SMP->True,
		TransversePolarizationVectors->{k},
		ExcludeParticles->{state\[Gamma]}
	]];
	m1=FeynArts`TheMass[inStates[[1]]];
	m2=FeynArts`TheMass[inStates[[2]]];
	m3=FeynArts`TheMass[outStates[[1]]];
	m4=FeynArts`TheMass[outStates[[2]]];
	(* set masses *)
	FeynCalc`SP[p1,p1]=m1^2;
	FeynCalc`SP[p2,p2]=m2^2;
	FeynCalc`SP[p3,p3]=m3^2;
	FeynCalc`SP[p4,p4]=m4^2;
	FeynCalc`SP[k,k]=0;
	
	(* strip off polarization vectors *)
	ampIS=ReplaceAll[ampIS,{
		FeynCalc`DiracGamma[FeynCalc`Momentum[FeynCalc`Polarization[P,-I]]]:>FeynCalc`GA[\[Mu]],
		FeynCalc`Pair[a__FeynCalc`Momentum,FeynCalc`Momentum[FeynCalc`Polarization[P, I]]]:>FeynCalc`Pair[a,FeynCalc`LorentzIndex[\[Mu]]],
		FeynCalc`Pair[FeynCalc`Momentum[FeynCalc`Polarization[P, I]],a__FeynCalc`Momentum]:>FeynCalc`Pair[a,FeynCalc`LorentzIndex[\[Mu]]]}
		];
	(* square amplitude and replace mu with new Lorentz index*)
	L\[Mu]\[Nu]=ampIS*ReplaceAll[FeynCalc`ComplexConjugate[ampIS],\[Mu]->\[Nu]];
	(* average over initial state fermions spins *)
	L\[Mu]\[Nu]=ReplaceAll[FeynCalc`FermionSpinSum[L\[Mu]\[Nu],FeynCalc`ExtraFactor->1/4],FeynCalc`DiracTrace->FeynCalc`TR];
	(* Compute L using L\[Mu]\[Nu] = (P^\[Mu]P^\[Nu]-g^\[Mu]\[Nu]P^2)L => L = -g^\[Mu]\[Nu]Subscript[L, \[Mu]\[Nu]]/3Q^2 *)
	L=Simplify[FeynCalc`Contract[L\[Mu]\[Nu]*FeynCalc`MT[\[Mu],\[Nu]]]]/(3Q^2);
	
	(* Final state Integrate[M(V\[Rule]finalstate)^2/4,phasespace] *)
	
	(* strip off polarization vectors *)
	ampFS=ReplaceAll[ampFS,{
		FeynCalc`DiracGamma[FeynCalc`Momentum[FeynCalc`Polarization[P,I]]]:>FeynCalc`GA[\[Mu]],
		FeynCalc`Pair[a__FeynCalc`Momentum,FeynCalc`Momentum[FeynCalc`Polarization[P, I]]]:>FeynCalc`Pair[a,FeynCalc`LorentzIndex[\[Mu]]],
		FeynCalc`Pair[FeynCalc`Momentum[FeynCalc`Polarization[P, I]],a__FeynCalc`Momentum]:>FeynCalc`Pair[a,FeynCalc`LorentzIndex[\[Mu]]]}
		];
	(* square amplitude and replace mu with new Lorentz index*)
	X\[Mu]\[Nu]=ampFS*ReplaceAll[FeynCalc`ComplexConjugate[ampFS],\[Mu]->\[Nu]];
	(* sum over final state fermions spins *)
	X\[Mu]\[Nu]=ReplaceAll[FeynCalc`FermionSpinSum[X\[Mu]\[Nu]],FeynCalc`DiracTrace->FeynCalc`TR];
	(* sum over final state photon polarizations *)
	X\[Mu]\[Nu]=FeynCalc`DoPolarizationSums[X\[Mu]\[Nu],k,0];
	(* Compute X using X\[Mu]\[Nu] = (P^\[Mu]P^\[Nu]-g^\[Mu]\[Nu]P^2)X => X = -g^\[Mu]\[Nu]Subscript[X, \[Mu]\[Nu]]/3Q^2 *)
	X=Simplify[FeynCalc`Contract[X\[Mu]\[Nu]*FeynCalc`MT[\[Mu],\[Nu]]]]/(3Q^2);
	
	(* Compute M^2 = L\[Mu]\[Nu]X\[Mu]\[Nu] = 3Q^2 X L and add vector propagator *)
	msqrd=X*L (3*Q^4)/ ((Q^2-Global`mv^2)^2+(Global`mv*Global`widthv)^2);
	
	(* Set the Mandelstam variables *)
	FeynCalc`SP[p3,p4]=(s-m3^2-m4^2)/2; (* s = (P - k)^2 *)
	FeynCalc`SP[k, p4]=(t-m4^2)/2; (* t = (P - p3)^2 *)
	FeynCalc`SP[k, p3]=(u-m3^2)/2; (* u = (P - p4)^2 *)
	FeynCalc`SP[p1,p2]=(Q^2-m1^2-m2^2)/2;

	msqrd=ReplaceAll[msqrd,{u->Q^2+m3^2+m4^2-s-t}];

	tbounds=ComputeMandelstamTBound[Q, 0, m3, m4, s, t];
	E1=(Q^2+m1^2-m2^2)/(2Q);
	tintegral=Integrate[msqrd,t,GenerateConditions->False];
	dndE=Simplify[((tintegral/.tbounds[[1]])-(tintegral/.tbounds[[2]])),Assumptions->{Q>0,Q>m1+m2,Q>m3+m4}];
	preFactor=(2*Q)/(4*Q Sqrt[E1^2-m1^2]) 1/(16*Q^2 (2\[Pi])^3)*1/\[Sigma]0;
	Simplify[
		ReplaceAll[
			(-1)preFactor*dndE,
			{s->Q^2(1 - 2 Global`E\[Gamma] / Q),Global`qe->Sqrt[4\[Pi] Global`\[Alpha]em]}],
		Assumptions->{Q>0, Q-(m3+m4)^2/Q > 2 * Global`E\[Gamma] > 0,m1>0,m2>0,m3>0,m4>0}]
];


HazmaComputeDNDE[inStates_,outStates_,Q_,OptionsPattern[]]:=Module[{newOutStates,msqrd,\[Sigma]0,p1,p2,p3,p4,k,m1,m2,m3,m4,s,t,u,E1,E3,tbounds,tintegral,dndE,preFactor},
	(* TODO: Add options *)
	If[$HazmaModel==="scalar",Return[ScalarMediatorComputeDNDE[inStates,outStates,Q]],Continue];
	(* TODO: Add options *)
	If[$HazmaModel==="vector",Return[VectorMediatorComputeDNDE[inStates,outStates,Q]],Continue];
	(* TODO: Throw error - can't do axial vector yet. *)
	Throw["invalid model"]
];


End[]


EndPackage[]
