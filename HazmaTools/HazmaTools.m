(* ::Package:: *)

(* ::Title:: *)
(*HazmaTools*)


(* ::Author:: *)
(*Authors: Adam Coogan[a] and Logan Morrison[b,c]*)


(* ::Affiliation:: *)
(*[a] GRAPPA, Institute of Physics, University of Amsterdam, 1098 XH Amsterdam, The Nether- lands*)
(*[b]Department of Physics, 1156 High St., University of California Santa Cruz, Santa Cruz, CA 95064, USA*)
(*[c]Santa Cruz Institute for Particle Physics, 1156 High St., Santa Cruz, CA 95064, USA*)


(* ::Section:: *)
(*Package Setup*)


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


Lepton::usage="FeynArts state for generic down-type lepton";
AntiLepton::usage="FeynArts state for generic anti-down-type lepton";
Electron::usage="FeynArts state for electron";
AntiElectron::usage="FeynArts state for anti-electron";
Muon::usage="FeynArts state for muon";
AntiMuon::usage="FeynArts state for anti-muon";
Tau::usage="FeynArts state for tau";
AntiTau::usage="FeynArts state for anti-tau";
Neutrino::usage="FeynArts state for generic left-handed neutrino";
AntiNeutrino::usage="FeynArts state for generic anit left-handed neutrino";
ElectronNeutrino::usage="FeynArts state for electron-type neutrino";
AntiElectronNeutrino::usage="FeynArts state for anti electron-type neutrino";
MuonNeutrino::usage="FeynArts state for muon-type neutrino";
AntiMuonNeutrino::usage="FeynArts state for anti muon-type neutrino";
TauNeutrino::usage="FeynArts state for tau-type neutrino";
AntiTauNeutrino::usage="FeynArts state for anti tau-type neutrino";
UpTypeQuark::usage="FeynArts state for generic up-type quark";
AntiUpTypeQuark::usage="FeynArts state for generic anti up-type quark";
UpQuark::usage="FeynArts state for up-quark";
AntiUpQuark::usage="FeynArts state for anti up-quark";
CharmQuark::usage="FeynArts state for charm-quark";
AntiCharmQuark::usage="FeynArts state for anti charm-quark";
TopQuark::usage="FeynArts state for top-quark";
AntiTopQuark::usage="FeynArts state for anti top-quark";
DownTypeQuark::usage="FeynArts state for generic down-type quark";
AntiDownTypeQuark::usage="FeynArts state for generic anti down-type quark";
DownQuark::usage="FeynArts state for down-quark";
AntiDownQuark::usage="FeynArts state for anti down-quark";
StrangeQuark::usage="FeynArts state for strange-quark";
AntiStrangeQuark::usage="FeynArts state for anti strange-quark";
BottomQuark::usage="FeynArts state for bottom-quark";
AntiBottomQuark::usage="FeynArts state for anti botoom-quark";
Photon::usage="FeynArts state for photon";
ZBoson::usage="FeynArts state for Z";
WBosonMinus::usage="FeynArts state for W-";
WBosonPlus::usage="FeynArts state for W+";
Higgs::usage="FeynArts state for Higgs";
ZGoldstone::usage="FeynArts state for Golstone associated with Z";
WMinusGoldstone::usage="FeynArts state for Golstone associated with W-";
WPlusGoldstone::usage="FeynArts state for Golstone associated with W+";
PhotonGhost::usage="FeynArts state for ghost associated with photon";
AntiPhotonGhost::usage="FeynArts state for anti-ghost associated with photon";
ZGhost::usage="FeynArts state for ghost associated with Z";
AntiZGhost::usage="FeynArts state for anti-ghost associated with Z";
WMinusGhost::usage="FeynArts state for ghost associated with W-";
AntiWMinusGhost::usage="FeynArts state for anti-ghost associated with W-";
WPlusGhost::usage="FeynArts state for ghost associated with W+";
AntiWPlusGhost::usage="FeynArts state for anti-ghost associated with W+";


HazmaGenerateDiagrams::usage="HazmaGenerateDiagrams[inStates,outStates] Generate `FeynArts` Feynmann diagrams for the process `inStates`->`outStates`.";


HazmaSubstitutionUnstableMediatorPropagator::usage="HazmaSubstitutionUnstableMediatorPropagator[exp] Substitue the mediator propagator for its unstable form in `exp`";


HazmaCreateAmplitude::usage="HazmaCreateAmplitude[inStates,outStates] creates an amplitude for `inStates`->`outStates` in the Scalar-Mediator model.";


HazmaCreateAmplitudeSquared::usage="HazmaCreateAmplitudeSquared[inStates,outStates] creates the squared spin-averaged amplitude for `inStates`->`outStates` in the Scalar-Mediator model.";


HazmaComputeCrossSection22::usage="HazmaComputeCrossSection22[inStates,outStates,Q] computes 2->2 cross section forprocess `inStates`->`outStates` given a center of mass energy `Q` in the Scalar-Mediator model.";


HazmaComputeDNDE::usage="HazmaComputeDNDE[inStates,outStates,Q] computes \[Gamma]-ray spectrum for `inStates`->`outStates` + \[Gamma] given a center of mass energy `Q` in the Scalar-Mediator model. Returns an expression in terms of ";


HazmaComputeWidth12::usage="HazmaComputeWidth12[inState,outStates] computes the 1->2 decay width for inState->outStates."


HazmaComputeWidth13::usage="HazmaComputeWidth13[inState,outStates] computes the 1->3 decay width for inState->outStates."


Begin["`Private`"]


Print["MeV-DM - Scalar Mediator: version 0.2.1"];
Print["by Adam Coogan and Logan A. Morrison"];


$HazmaModel="scalar";


(* ::Section::Closed:: *)
(*Utilities*)


(* ::Subsection::Closed:: *)
(*Misc.*)


ComputeMandelstamTBound[M_, m1_, m2_, m3_, s_, t_]:=Module[{E1,E2},
	E1=(M^2+m1^2-s)/(2*M);
	E2=(M^2+m2^2-t)/(2*M);
	Solve[M^2+2*E1*E2+m1^2+m2^2-m3^2-2*M*(E1+E2)==2Sqrt[E1^2-m1^2] Sqrt[E2^2-m2^2],t]
];


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


(* SM States from "SM" model for testing *)
(* down-type leptons *)
Lepton=FeynArts`F[2];
AntiLepton=FeynArts`F[2];
Electron=FeynArts`F[2,{1}];
AntiElectron=-FeynArts`F[2,{1}];
Muon=FeynArts`F[2,{2}];
AntiMuon=-FeynArts`F[2,{2}];
Tau=FeynArts`F[2,{3}];
AntiTau=-FeynArts`F[2,{3}];
(* neutinos *)
Neutrino=FeynArts`F[1];
AntiNeutrino=FeynArts`F[1];
ElectronNeutrino=FeynArts`F[1,{1}];
AntiElectronNeutrino=-FeynArts`F[1,{1}];
MuonNeutrino=FeynArts`F[1,{2}];
AntiMuonNeutrino=-FeynArts`F[1,{2}];
TauNeutrino=FeynArts`F[1,{3}];
AntiTauNeutrino=-FeynArts`F[1,{3}];
(* up-type quarks *)
UpTypeQuark=FeynArts`F[3];
AntiUpTypeQuark=FeynArts`F[3];
UpQuark=FeynArts`F[3,{1}];
AntiUpQuark=-FeynArts`F[3,{1}];
CharmQuark=FeynArts`F[3,{2}];
AntiCharmQuark=-FeynArts`F[3,{2}];
TopQuark=FeynArts`F[3,{3}];
AntiTopQuark=-FeynArts`F[3,{3}];
(* down-type quarks *)
DownTypeQuark=FeynArts`F[4];
AntiDownTypeQuark=FeynArts`F[4];
DownQuark=FeynArts`F[4,{1}];
AntiDownQuark=-FeynArts`F[4,{1}];
StrangeQuark=FeynArts`F[4,{2}];
AntiStrangeQuark=-FeynArts`F[4,{2}];
BottomQuark=FeynArts`F[4,{3}];
AntiBottomQuark=-FeynArts`F[4,{3}];
(* vector bosons *)
Photon=FeynArts`V[1];
ZBoson=FeynArts`V[2];
WBosonMinus=FeynArts`V[3];
WBosonPlus=-FeynArts`V[3];
(* Higgs *)
Higgs=FeynArts`S[1];
(* Golstones *)
ZGoldstone=FeynArts`S[2];
WMinusGoldstone=FeynArts`S[3];
WPlusGoldstone=-FeynArts`S[3];
(* ghosts *)
PhotonGhost=FeynArts`U[1];
AntiPhotonGhost=FeynArts`U[1];
ZGhost=FeynArts`U[2];
AntiZGhost=-FeynArts`U[2];
WMinusGhost=FeynArts`U[3];
AntiWMinusGhost=-FeynArts`U[3];
WPlusGhost=FeynArts`U[4];
AntiWPlusGhost=-FeynArts`U[4];


(* determine if state is a fermion *)
isFermion[state_]:=Module[{head},
	(* if state looks like -F[1], need to strip off the `Times` head *)
	head=If[Head[state]===Times, Head[ReplaceAll[state,Times[a_,b__]:>b]], Head[state]];
	head===FeynArts`F \[Or] head===FeynArts`U
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
	If[$HazmaModel==="SM", Return["SM"]]; (*for testing*)
	If[$HazmaModel==="SMQCD", Return["SMQCD"]]; (*for testing*)
	Throw[$Failed, Abort::InvalidHazmaModel]
]


GetGenericModelName[]:=Module[{modelName},
	If[$HazmaModel==="scalar",Return[{"Lorentz","EFT_MeV_DM_scalar"<>"/"<>"EFT_MeV_DM_scalar"}]];
	If[$HazmaModel==="vector",Return[{"Lorentz","EFT_MeV_DM_vector_no_contact"<>"/"<>"EFT_MeV_DM_vector_no_contact"}]];
	If[$HazmaModel==="SM", Return[OptionValue[FeynArts`InsertFields,FeynArts`GenericModel]]]; (*for testing*)
	If[$HazmaModel==="SMQCD", Return[OptionValue[FeynArts`InsertFields,FeynArts`GenericModel]]]; (*for testing*)
	Throw[$Failed, Abort::InvalidHazmaModel]
]


Options[HazmaGenerateDiagrams]={
	FeynArts`LoopNumber->0,
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`ExcludeParticles->{},
	FeynArts`ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	FeynArts`InsertionLevel->{FeynArts`Particles},
	FeynArts`Paint->False
};

HazmaGenerateDiagrams[inStates_, outStates_,OptionsPattern[]]:=Block[{tops,diags},
	tops = FeynArts`CreateTopologies[
		OptionValue[FeynArts`LoopNumber], 
		Length[inStates]-> Length[outStates],
		FeynArts`ExcludeTopologies-> OptionValue[FeynArts`ExcludeTopologies],
		FeynArts`Adjacencies-> OptionValue[FeynArts`Adjacencies]
	];
	diags=FeynArts`InsertFields[
		tops,
		inStates -> outStates,
		FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
		FeynArts`Model -> GetModelName[],
		FeynArts`GenericModel -> GetGenericModelName[],
		FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles]
	];
	If[OptionValue[FeynArts`Paint],FeynArts`Paint[diags],FeynArts`SheetHeader->None];
	diags
]


(* ::Section::Closed:: *)
(*Amplitude Creation*)


(* ::Subsection:: *)
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
	If[$HazmaModel==="SM", Return[exp]];
	If[$HazmaModel==="SMQCD", Return[exp]];
	Throw[$Failed, Abort::InvalidHazmaModel]
]


Clear[HazmaCreateAmplitude];

Options[HazmaCreateAmplitude]={
	FeynArts`LoopNumber->0,
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`ExcludeParticles->{},
	FeynArts`ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	FeynArts`InsertionLevel->{FeynArts`Particles},
	FeynArts`Paint->False,
	FeynCalc`ChangeDimension->4,
	FeynCalc`FinalSubstitutions->{Global`M$FACouplings},
	FeynCalc`IncomingMomenta->{},
	List->False,
	FeynCalc`LoopMomenta->{},
	FeynCalc`LorentzIndexNames->{Global`\[Mu],Global`\[Nu],Global`\[Alpha],Global`\[Beta],Global`\[Rho],Global`\[Sigma]},
	FeynCalc`OutgoingMomenta->{},
	FeynCalc`SMP->True,
	FeynCalc`TransversePolarizationVectors->{},
	FeynCalc`UndoChiralSplittings->True,
	FeynArts`GaugeRules-> _FeynArts`FAGaugeXi->1,
	FeynArts`PreFactor->-I (2 \[Pi])^(-4 * FeynArts`LoopNumber),
	FeynArts`Truncated->False
};

HazmaCreateAmplitude[inStates_,outStates_,OptionsPattern[]] := Module[{tops, ampFA,ampFC, diagrams},

FeynCalc`ClearScalarProducts[];
(* create the topologies based on input states *)
diagrams=HazmaGenerateDiagrams[inStates,outStates,
	FeynArts`LoopNumber->OptionValue[FeynArts`LoopNumber],
	FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
	FeynArts`ExcludeParticles->OptionValue[FeynArts`ExcludeParticles],
	FeynArts`ExcludeTopologies->OptionValue[FeynArts`ExcludeTopologies],
	FeynArts`InsertionLevel->OptionValue[FeynArts`InsertionLevel],
	FeynArts`Paint->OptionValue[FeynArts`Paint]
];
	
(*Compute the amplitude *)
ampFA=FeynArts`CreateFeynAmp[diagrams, 
	FeynArts`GaugeRules->OptionValue[FeynArts`GaugeRules],
	FeynArts`PreFactor->OptionValue[FeynArts`PreFactor],
	FeynArts`Truncated->OptionValue[FeynArts`Truncated]
];

ampFC=FeynCalc`FCFAConvert[
	ampFA,
	FeynCalc`IncomingMomenta-> OptionValue[FeynCalc`IncomingMomenta],
	FeynCalc`OutgoingMomenta-> OptionValue[FeynCalc`OutgoingMomenta],
	FeynCalc`LorentzIndexNames-> OptionValue[FeynCalc`LorentzIndexNames],
	FeynCalc`LoopMomenta->OptionValue[FeynCalc`LoopMomenta], 
	FeynCalc`UndoChiralSplittings->OptionValue[FeynCalc`UndoChiralSplittings],
	FeynCalc`SMP-> OptionValue[FeynCalc`SMP],
	FeynCalc`ChangeDimension-> OptionValue[FeynCalc`ChangeDimension],
	FeynCalc`FinalSubstitutions-> OptionValue[FeynCalc`FinalSubstitutions],
	FeynCalc`TransversePolarizationVectors->OptionValue[FeynCalc`TransversePolarizationVectors]
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


(* ::Subsection::Closed:: *)
(*Computing Squared Amplitudes*)


Options[HazmaCreateAmplitudeSquared]={
	FeynArts`LoopNumber->0,
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`ExcludeParticles->{},
	FeynArts`ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	FeynArts`InsertionLevel->{FeynArts`Particles},
	FeynArts`Paint->False,
	FeynArts`GaugeRules->_FeynArts`FAGaugeXi->1,
	FeynArts`PreFactor->-I (2 \[Pi])^(-4 FeynArts`LoopNumber),
	FeynArts`Truncated->False,
	FeynCalc`ChangeDimension->4,
	FeynCalc`FinalSubstitutions->{Global`M$FACouplings},
	FeynCalc`IncomingMomenta->{},
	List->False,
	FeynCalc`LoopMomenta->{},
	FeynCalc`LorentzIndexNames->{Global`\[Mu],Global`\[Nu],Global`\[Alpha],Global`\[Beta],Global`\[Rho],Global`\[Sigma]},
	FeynCalc`OutgoingMomenta->{},
	FeynCalc`SMP->True,
	FeynCalc`TransversePolarizationVectors->{},
	FeynCalc`UndoChiralSplittings->True
};

HazmaCreateAmplitudeSquared[inStates_, outStates_,OptionsPattern[]] := Module[{msqrd, amp,dim,momentum,inMomenta,outMomenta},

(* Create the amplitude *)
amp=HazmaCreateAmplitude[
	inStates, 
	outStates,
	FeynArts`LoopNumber->OptionValue[FeynArts`LoopNumber],
	FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
	FeynArts`ExcludeParticles->OptionValue[FeynArts`ExcludeParticles],
	FeynArts`ExcludeTopologies->OptionValue[FeynArts`ExcludeTopologies],
	FeynArts`InsertionLevel->OptionValue[FeynArts`InsertionLevel],
	FeynArts`Paint->OptionValue[FeynArts`Paint],
	FeynArts`GaugeRules->OptionValue[FeynArts`GaugeRules],
	FeynArts`PreFactor->OptionValue[FeynArts`PreFactor],
	FeynArts`Truncated->OptionValue[FeynArts`Truncated],
	FeynCalc`ChangeDimension->OptionValue[FeynCalc`ChangeDimension],
	FeynCalc`FinalSubstitutions->OptionValue[FeynCalc`FinalSubstitutions],
	FeynCalc`IncomingMomenta->OptionValue[FeynCalc`IncomingMomenta],
	List->OptionValue[List],
	FeynCalc`LoopMomenta->OptionValue[FeynCalc`LoopMomenta],
	FeynCalc`LorentzIndexNames->OptionValue[FeynCalc`LorentzIndexNames],
	FeynCalc`OutgoingMomenta->OptionValue[FeynCalc`OutgoingMomenta],
	FeynCalc`SMP->OptionValue[FeynCalc`SMP],
	FeynCalc`TransversePolarizationVectors->OptionValue[FeynCalc`TransversePolarizationVectors],
	FeynCalc`UndoChiralSplittings->OptionValue[FeynCalc`UndoChiralSplittings]
	];

FeynCalc`ClearScalarProducts[];
(* determine incoming/outgoing momenta *)
inMomenta=OptionValue[FeynCalc`IncomingMomenta];
outMomenta=OptionValue[FeynCalc`OutgoingMomenta];

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


(* ::Section::Closed:: *)
(*Cross Sections*)


(* ::Subsection::Closed:: *)
(*Compute 2->2 Cross Sections*)


Options[HazmaComputeCrossSection22]={
	FeynArts`LoopNumber->0,
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`ExcludeParticles->{},
	FeynArts`ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	FeynArts`InsertionLevel->{FeynArts`Particles},
	FeynArts`GaugeRules->_FeynArts`FAGaugeXi->1,
	FeynArts`Paint->False,
	FeynCalc`ChangeDimension->4,
	FeynCalc`FinalSubstitutions->{Global`M$FACouplings},
	FeynCalc`LoopMomenta->{},
	FeynCalc`SMP->True
};

HazmaComputeCrossSection22[inStates_, outStates_,Q_,OptionsPattern[]] := Module[{msqrd,p1,p2,p3,p4,m1,m2,m3,m4,s,t,u,E1,E3,pi,pf,preFactor,tmin,tmax},
	ClearScalarProducts[];
	(* Create the amplitude *)
	msqrd=HazmaCreateAmplitudeSquared[
		inStates, 
		outStates,
		FeynArts`LoopNumber->OptionValue[FeynArts`LoopNumber],
		FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
		FeynArts`ExcludeParticles->OptionValue[FeynArts`ExcludeParticles],
		FeynArts`ExcludeTopologies->OptionValue[FeynArts`ExcludeTopologies],
		FeynArts`InsertionLevel->OptionValue[FeynArts`InsertionLevel],
		FeynArts`Paint->OptionValue[FeynArts`Paint],
		FeynArts`GaugeRules->OptionValue[FeynArts`GaugeRules],
		FeynCalc`ChangeDimension->OptionValue[FeynCalc`ChangeDimension],
		FeynCalc`FinalSubstitutions->OptionValue[FeynCalc`FinalSubstitutions],
		FeynCalc`IncomingMomenta->{p1,p2},
		FeynCalc`LoopMomenta->OptionValue[FeynCalc`LoopMomenta],
		FeynCalc`OutgoingMomenta->{p3,p4},
		FeynCalc`SMP->OptionValue[FeynCalc`SMP]
	];
	
	m1=FeynArts`TheMass[inStates[[1]]];
	m2=FeynArts`TheMass[inStates[[2]]];
	m3=FeynArts`TheMass[outStates[[1]]];
	m4=FeynArts`TheMass[outStates[[2]]];

	FeynCalc`SetMandelstam[s,t,u,p1,p2,-p3,-p4,m1,m2,m3,m4];	
	msqrd=Simplify[msqrd];
	msqrd=ReplaceAll[msqrd,{u->m1^2+m2^2+m3^2+m4^2-s-t}];

	E1=(s+m1^2-m2^2)/(2*Sqrt[s]);
	E3=(s+m3^2-m4^2)/(2*Sqrt[s]);
	pi=Sqrt[E1^2-m1^2]; (* mag. of initial state 3-momentum *)
	pf=Sqrt[E3^2-m3^2]; (* mag. of final state 3-momentum *)
	(*msqrd=Simplify[ReplaceAll[msqrd,{t->m1^2+m3^2-2(E1*E3-z p1*p3)}]];*)
	
	preFactor=Simplify[(2*\[Pi])/(64*\[Pi]^2*Q^2)*(pf/pi)*(1/(2*pi*pf))];
	(* If final state particles are identical, divide by 2 *)
	If[outStates[[1]]===outStates[[2]],preFactor=preFactor/2,Continue];
	
	tmax=Simplify[m1^2+m3^2-2E1*E3+2*pi*pf];
	tmin=Simplify[m1^2+m3^2-2E1*E3-2*pi*pf];

	ReplaceAll[preFactor*Integrate[msqrd,{t,tmin,tmax},GenerateConditions->False],{s->Q^2}]

];


(* ::Section::Closed:: *)
(*Widths*)


(* ::Subsection:: *)
(*Compute Width 1->2*)


Options[HazmaComputeWidth12]={
	FeynArts`LoopNumber->0,
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`ExcludeParticles->{},
	FeynArts`ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	FeynArts`InsertionLevel->{FeynArts`Particles},
	FeynArts`Paint->False,
	FeynArts`GaugeRules->_FeynArts`FAGaugeXi->1,
	FeynCalc`ChangeDimension->4,
	FeynCalc`FinalSubstitutions->{Global`M$FACouplings},
	FeynCalc`LoopMomenta->{},
	FeynCalc`SMP->True
};

HazmaComputeWidth12[inState_, outStates_,OptionsPattern[]] := Module[{msqrd,P,p1,p2,M,m1,m2,pf,preFactor,E1,E2},
	ClearScalarProducts[];
	If[Length[inState]!=1,Throw[$Failed,HazmaComputeWidth12::InvalidInState],Continue];
	If[Length[outStates]!=2,Throw[$Failed,HazmaComputeWidth12::InvalidOutState],Continue];
	(* Create the amplitude *)
	msqrd=HazmaCreateAmplitudeSquared[
		{inState}, 
		outStates,
		FeynArts`LoopNumber->OptionValue[FeynArts`LoopNumber],
		FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
		FeynArts`ExcludeParticles->OptionValue[FeynArts`ExcludeParticles],
		FeynArts`ExcludeTopologies->OptionValue[FeynArts`ExcludeTopologies],
		FeynArts`InsertionLevel->OptionValue[FeynArts`InsertionLevel],
		FeynArts`Paint->OptionValue[FeynArts`Paint],
		FeynArts`GaugeRules->OptionValue[FeynArts`GaugeRules],
		FeynCalc`ChangeDimension->OptionValue[FeynCalc`ChangeDimension],
		FeynCalc`FinalSubstitutions->OptionValue[FeynCalc`FinalSubstitutions],
		FeynCalc`IncomingMomenta->{P},
		FeynCalc`LoopMomenta->OptionValue[FeynCalc`LoopMomenta],
		FeynCalc`OutgoingMomenta->{p1,p2},
		FeynCalc`SMP->OptionValue[FeynCalc`SMP]
	];
	
	M=FeynArts`TheMass[inState];
	m1=FeynArts`TheMass[outStates[[1]]];
	m2=FeynArts`TheMass[outStates[[2]]];
	
	E1=(M^2+m1^2-m2^2)/(2*M);
	E2=(M^2-m1^2+m2^2)/(2*M);
	
	FeynCalc`SP[P,P]=M^2;
	FeynCalc`SP[p1,p1]=m1^2;	
	FeynCalc`SP[p2,p2]=m2^2;
	FeynCalc`SP[P,p1]=M*E1;	
	FeynCalc`SP[P,p2]=M*E2;
	FeynCalc`SP[p1,p2]=(M^2-m1^2-m2^2)/2;
	
	msqrd=Simplify[msqrd];

	pf=Sqrt[E1^2-m1^2]; (* mag. of final state 3-momentum *)
	
	preFactor=Simplify[(4*\[Pi])/(2*M)*pf/(16*Pi^2*M)];
	(* If final state particles are identical, divide by 2 *)
	preFactor*=If[outStates[[1]]===outStates[[2]],1/2, 1];
	preFactor*msqrd

];


(* ::Subsection:: *)
(*Compute Width 1->3*)


Options[HazmaComputeWidth13]={
	FeynArts`LoopNumber->0,
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`ExcludeParticles->{},
	FeynArts`ExcludeTopologies->{FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
	FeynArts`InsertionLevel->{FeynArts`Particles},
	FeynArts`Paint->False,
	FeynArts`GaugeRules->_FeynArts`FAGaugeXi->1,
	FeynCalc`ChangeDimension->4,
	FeynCalc`FinalSubstitutions->{Global`M$FACouplings},
	FeynCalc`LoopMomenta->{},
	FeynCalc`SMP->True
};

HazmaComputeWidth13[inState_, outStates_,OptionsPattern[]] := Module[{msqrd,P,p1,p2,p3,M,m1,m2,m3,pf,preFactor,tbounds,tintegral,sintegral},
	ClearScalarProducts[];
	If[Length[inState]!=1,Throw[$Failed,HazmaComputeWidth13::InvalidInState],Continue];
	If[Length[outStates]!=3,Throw[$Failed,HazmaComputeWidth13::InvalidOutState],Continue];
	(* Create the amplitude *)
	msqrd=HazmaCreateAmplitudeSquared[
		{inState}, 
		outStates,
		FeynArts`LoopNumber->OptionValue[FeynArts`LoopNumber],
		FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
		FeynArts`ExcludeParticles->OptionValue[FeynArts`ExcludeParticles],
		FeynArts`ExcludeTopologies->OptionValue[FeynArts`ExcludeTopologies],
		FeynArts`InsertionLevel->OptionValue[FeynArts`InsertionLevel],
		FeynArts`Paint->OptionValue[FeynArts`Paint],
		FeynArts`GaugeRules->OptionValue[FeynArts`GaugeRules],
		FeynCalc`ChangeDimension->OptionValue[FeynCalc`ChangeDimension],
		FeynCalc`FinalSubstitutions->OptionValue[FeynCalc`FinalSubstitutions],
		FeynCalc`IncomingMomenta->{P},
		FeynCalc`LoopMomenta->OptionValue[FeynCalc`LoopMomenta],
		FeynCalc`OutgoingMomenta->{p1,p2,p3},
		FeynCalc`SMP->OptionValue[FeynCalc`SMP]
	];
	
	M=FeynArts`TheMass[inState];
	m1=FeynArts`TheMass[outStates[[1]]];
	m2=FeynArts`TheMass[outStates[[2]]];
	m3=FeynArts`TheMass[outStates[[3]]];
	
	FeynCalc`SetMandelstam[s,t,u,P,-p1,-p2,-p3,M,m1,m2,m3];
	msqrd=Simplify[msqrd];
	msqrd=ReplaceAll[msqrd,u->{M^2+m1^2+m2^2+m3^2-s-t}];
	
	tbounds=ComputeMandelstamTBound[M,m1,m2,m3,s,t];
	
	(* Integrate over t *)
	tintegral=Integrate[msqrd,t,GenerateConditions->False];
	tintegral=ReplaceAll[tintegral,tbounds[[2]]]-ReplaceAll[tintegral,tbounds[[1]]];
	(* Integrate over s *)
	sintegral=Integrate[tintegral,s,GenerateConditions->False];
	sintegral=ReplaceAll[sintegral,s->(M-m1)^2]-ReplaceAll[sintegral,s->(m2+m3)^2];
	
	preFactor=Simplify[1/(2*M)*1/(16*M^2*(2*Pi)^3)];
	(* If final state particles are identical, divide by 2 *)
	preFactor/=Factorial[4-CountDistinct[outStates]];
	preFactor*sintegral

];


(* ::Section::Closed:: *)
(*dNdE*)


(* ::Subsection:: *)
(*Computing dNdE*)


Options[ScalarMediatorComputeDNDE]={
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`Paint->False,
	FeynCalc`FinalSubstitutions->{Global`M$FACouplings}
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
		FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
		FeynArts`Paint->OptionValue[FeynArts`Paint],
		(* ignore diagrams with internal photon. They are higher order in \[Alpha] *)
		FeynArts`ExcludeParticles->{state\[Gamma]},
		FeynCalc`ChangeDimension->4,
		FeynCalc`FinalSubstitutions->Global`M$FACouplings,
		FeynCalc`IncomingMomenta->{p1,p2},
		List->False,
		FeynCalc`LorentzIndexNames->{\[Mu],\[Nu],\[Alpha],\[Beta]},
		FeynCalc`OutgoingMomenta->{p3,p4,k},
		FeynCalc`SMP->True,
		FeynCalc`TransversePolarizationVectors->{k}
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
	FeynArts`Adjacencies->{3,4,5},
	FeynArts`Paint->False,
	FeynCalc`FinalSubstitutions->{Global`M$FACouplings}
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
		FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
		FeynArts`Paint->OptionValue[FeynArts`Paint],
		FeynArts`ExcludeParticles->{state\[Gamma]},
		FeynCalc`ChangeDimension->4,
		FeynCalc`FinalSubstitutions->Global`M$FACouplings,
		FeynCalc`IncomingMomenta->{p1,p2},
		List->False,
		FeynCalc`LorentzIndexNames->{},
		FeynCalc`OutgoingMomenta->{P},
		FeynCalc`SMP->True,
		FeynCalc`TransversePolarizationVectors->{}
	]];
	(* Final state piece *)
	ampFS=FeynCalc`Contract[HazmaCreateAmplitude[
		{statev}, 
		newOutStates,
		FeynArts`Adjacencies->OptionValue[FeynArts`Adjacencies],
		FeynArts`Paint->OptionValue[FeynArts`Paint],
		FeynArts`ExcludeParticles->{state\[Gamma]},
		FeynCalc`ChangeDimension->4,
		FeynCalc`FinalSubstitutions->Global`M$FACouplings,
		FeynCalc`IncomingMomenta->{P},
		List->False,
		FeynCalc`LorentzIndexNames->{},
		FeynCalc`OutgoingMomenta->{p3,p4,k},
		FeynCalc`SMP->True,
		FeynCalc`TransversePolarizationVectors->{k}
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
