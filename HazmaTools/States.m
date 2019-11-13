(* Wolfram Language package *)

BeginPackage["HazmaTools`"]

NeutralPion::usage="FeynArts state for neutral pion";
ChargedPionM::usage="FeynArts state for negatively charged pion";
ChargedPionP::usage="FeynArts state for positively charged pion";
NeutralKaon::usage="FeynArts state for neutral kaon";
AntiNeutralKaon::usage="FeynArts state for anti neutral kaon";
ChargedKaonM::usage="FeynArts state for negatively charged kaon";
ChargedKaonP::usage="FeynArts state for positively charged kaon";
Eta::usage="FeynArts state for eta";
Rho::usage="FeynArts state for rho";
Omega::usage="FeynArts state for omega";
DarkMatter::usage="FeynArts state for dark matter";
AntiDarkMatter::usage="FeynArts state for anti dark matter";
ScalarMediator::usage="FeynArts state for scalar-mediator";
PseudoScalarMediator::usage="FeynArts state for pseudo-scalar-mediator";
VectorMediator::usage="FeynArts state for vector-mediator";


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
WBosonM::usage="FeynArts state for W-";
WBosonP::usage="FeynArts state for W+";
Gluon::usage="FeynArts state for gluons";
Higgs::usage="FeynArts state for Higgs";
GoldstoneZBoson::usage="FeynArts state for Golstone associated with Z";
GoldstoneWBosonM::usage="FeynArts state for Golstone associated with W-";
GoldstoneWBosonP::usage="FeynArts state for Golstone associated with W+";
GhostPhoton::usage="FeynArts state for ghost associated with photon";
AntiGhostPhoton::usage="FeynArts state for anti-ghost associated with photon";
GhostZBoson::usage="FeynArts state for ghost associated with Z";
AntiGhostZBoson::usage="FeynArts state for anti-ghost associated with Z";
GhostWBosonM::usage="FeynArts state for ghost associated with W-";
AntiGhostWBosonM::usage="FeynArts state for anti-ghost associated with W-";
GhostWBosonP::usage="FeynArts state for ghost associated with W+";
AntiGhostWBosonP::usage="FeynArts state for anti-ghost associated with W+";

StateInModelQ::usage = "StateInModelQ[state] Returns True if `state` is in \
current `$HazmaModel`"

AssertStateInModel::usage = "AssertStateInModel[state] Checks whether or not state \
is in the current $HazmaModel. Error is thrown if not."

isFermion::usage = "isFermion[state] Determine if a state is a fermion or not."
isVector::usage = "isVector[state] Determine if a state is a vector or not."

Begin["Private`"]

(* determine if state is a fermion *)
isFermion[state_]:=Module[{head},
	(* if state looks like -F[1], need to strip off the `Times` head *)
	head=If[Head[state]===Times, Head[ReplaceAll[state,Times[a_,b__]:>b]], Head[state]];
	head===FeynArts`F || head===FeynArts`U
]

(* determine if state is a vector *)
isVector[state_]:=Module[{head},
	(* if state looks like -V[1], need to strip off the `Times` head *)
	head=If[Head[state]===Times, Head[ReplaceAll[state,Times[a_,b__]:>b]], Head[state]];
	head===FeynArts`V
]	

StateInModelQ[NeutralPion] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[ChargedPionM] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[ChargedPionP] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[NeutralKaon] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiNeutralKaon] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[ChargedKaonM] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[ChargedKaonP] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[Eta] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];

StateInModelQ[Photon] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[Rho] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[Omega] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[Neutrino] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiNeutrino] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[Lepton] := MemberQ[{"scalar", "vector", "SM", "SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiLepton] := MemberQ[{"scalar", "vector", "SM", "SMQCD"}, HazmaTools`$HazmaModel];

StateInModelQ[DarkMatter] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiDarkMatter] := MemberQ[{"scalar", "vector"}, HazmaTools`$HazmaModel];

StateInModelQ[ScalarMediator] := MemberQ[{"scalar"}, HazmaTools`$HazmaModel];
StateInModelQ[PseudoScalarMediator] := MemberQ[{"pseudo-scalar"}, HazmaTools`$HazmaModel];
StateInModelQ[VectorMediator] := MemberQ[{"vector"}, HazmaTools`$HazmaModel];

StateInModelQ[Electron] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiElectron] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[Muon] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiMuon] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[Tau] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiTau] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];

StateInModelQ[ElectronNeutrino] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiElectronNeutrino] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[MuonNeutrino] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiMuonNeutrino] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[TauNeutrino] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiTauNeutrino] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];

StateInModelQ[UpTypeQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiUpTypeQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[UpQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiUpQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[CharmQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiCharmQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[TopQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiTopQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];

StateInModelQ[DownTypeQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiDownTypeQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[DownQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiDownQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[StrangeQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiStrangeQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[BottomQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiBottomQuark] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];

StateInModelQ[ZBoson] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[WBosonM] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[WBosonP] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[Gluon] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];

StateInModelQ[Higgs] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[GoldstoneZBoson] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[GoldstoneWBosonM] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[GoldstoneWBosonP] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];

StateInModelQ[GhostPhoton] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiGhostPhoton] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[GhostZBoson] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiGhostZBoson] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[GhostWBosonM] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiGhostWBosonM] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[GhostWBosonP] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];
StateInModelQ[AntiGhostWBosonP] := MemberQ[{"SM","SMQCD"}, HazmaTools`$HazmaModel];

AssertStateInModel::StateNotInModel = "The state `1 is not in the model `2";

AssertStateInModel[state_] := If[
	StateInModelQ[state],
	Continue,
	Message[AssertStateInModel::StateNotInModel, state, HazmaTools`$HazmaModel];
	Throw[$Failed]
];



(* SM Mesons *)
NeutralPion=FeynArts`S[1];
ChargedPionM=FeynArts`S[2];
ChargedPionP=-FeynArts`S[2];
NeutralKaon=FeynArts`S[3];
AntiNeutralKaon=-FeynArts`S[3];
ChargedKaonM=FeynArts`S[4];
ChargedKaonP=-FeynArts`S[4];
Eta=FeynArts`S[5];

(* SM vectors *)
Photon=FeynArts`V[1];
Rho=FeynArts`V[2];
Omega=FeynArts`V[1];

(* Leptons *)
Neutrino=FeynArts`F[1];
AntiNeutrino=-FeynArts`V[1];
Lepton=FeynArts`F[2];
AntiLepton=-FeynArts`F[2];

(* DM *)
DarkMatter=FeynArts`F[3];
AntiDarkMatter=-FeynArts`F[3];

(* mediators *)
ScalarMediator=FeynArts`S[7];
PseudoScalarMediator=FeynArts`S[7];
VectorMediator=FeynArts`V[4];

(* SM/SMQCD only states *)

(* down-type leptons *)
Electron=FeynArts`F[2,{1}];
AntiElectron=-FeynArts`F[2,{1}];
Muon=FeynArts`F[2,{2}];
AntiMuon=-FeynArts`F[2,{2}];
Tau=FeynArts`F[2,{3}];
AntiTau=FeynArts`F[2,{3}];

(* neutinos *)
ElectronNeutrino=FeynArts`F[1,{1}];
AntiElectronNeutrino=-FeynArts`F[1,{1}];
MuonNeutrino=FeynArts`F[1,{2}];
AntiMuonNeutrino=-FeynArts`F[1,{2}];
TauNeutrino=FeynArts`F[1,{3}];
AntiTauNeutrino=-FeynArts`F[1,{3}];

(* up-type quarks *)
UpTypeQuark=FeynArts`F[3];
AntiUpTypeQuark=-FeynArts`F[3];
UpQuark=FeynArts`F[3,{1}];
AntiUpQuark=-FeynArts`F[3,{1}];
CharmQuark=FeynArts`F[3,{2}];
AntiCharmQuark=-FeynArts`F[3,{2}];
TopQuark=FeynArts`F[3,{3}];
AntiTopQuark=-FeynArts`F[3,{3}];

(* down-type quarks *)
DownTypeQuark=FeynArts`F[4];
AntiDownTypeQuark=-FeynArts`F[4];
DownQuark=FeynArts`F[4,{1}];
AntiDownQuark=-FeynArts`F[4,{1}];
StrangeQuark=FeynArts`F[4,{2}];
AntiStrangeQuark=-FeynArts`F[4,{2}];
BottomQuark=FeynArts`F[4,{3}];
AntiBottomQuark=-FeynArts`F[4,{3}];

(* vector bosons *)
ZBoson=FeynArts`V[2];
WBosonM=FeynArts`V[3];
WBosonP=-FeynArts`V[3];
Gluon=-FeynArts`V[4];

(* Higgs and Goldstones *)
Higgs=FeynArts`S[1];
GoldstoneZBoson=FeynArts`S[2];
GoldstoneWBosonM=FeynArts`S[3];
GoldstoneWBosonP=-FeynArts`S[3];

(* ghosts *)
GhostPhoton=FeynArts`U[1];
AntiGhostPhoton=-FeynArts`U[1];
GhostZBoson=FeynArts`U[2];
AntiGhostZBoson=-FeynArts`U[2];
GhostWBosonM=FeynArts`U[3];
AntiGhostWBosonM=-FeynArts`U[3];
GhostWBosonP=FeynArts`U[4];
AntiGhostWBosonP=-FeynArts`U[5];

End[]


EndPackage[]