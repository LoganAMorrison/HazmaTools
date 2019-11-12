(* Mathematica Package *)

(* :Title: Amplitudes *)
(* :Context: HazmaTools` *)
(* :Author: Adam Coogan and Logan Morrison *)
(* :Date: 2019-11-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 Adam Coogan and Logan Morrison *)
(* :Keywords: Hazma, Quantum Field Theory, FeynCalc, FeynArts *)
(* :Discussion: File for computing amplitudes in the Hazma framework *)

BeginPackage["HazmaTools`"]

HazmaComputeAmplitude::usage = "HazmaComputeAmplitude[inState,outState] Compute the amplitude for \
inState->outState"
HazmaComputeAmplitudeSquared::usage = "HazmaComputeAmplitudeSquared[inState,outState] Compute the \
squared amplitude for inState->outState"

Begin["`Private`"]

Options[HazmaComputeAmplitude] := {
  FeynArts`LoopNumber -> 0,
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`ExcludeTopologies -> {FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`Paint -> False,
  FeynCalc`ChangeDimension -> 4,
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynCalc`IncomingMomenta -> {},
  List -> False,
  FeynCalc`LoopMomenta -> {},
  FeynCalc`LorentzIndexNames -> {Global`\[Mu], Global`\[Nu], Global`\[Alpha], Global`\[Beta], Global`\[Rho], Global`\[Sigma]},
  FeynCalc`OutgoingMomenta -> {},
  FeynCalc`SMP -> True,
  FeynCalc`TransversePolarizationVectors -> {},
  FeynCalc`UndoChiralSplittings -> True,
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynArts`PreFactor -> -I (2 \[Pi])^(-4 * FeynArts`LoopNumber),
  FeynArts`Truncated -> False
};

HazmaComputeAmplitude[inStates_, outStates_, OptionsPattern[]] := Module[{ampFA, ampFC, diagrams},

  FeynCalc`ClearScalarProducts[];
  (* create the topologies based on input states *)
  diagrams = HazmaGenerateDiagrams[inStates, outStates,
    FeynArts`LoopNumber -> OptionValue[FeynArts`LoopNumber],
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynArts`ExcludeTopologies -> OptionValue[FeynArts`ExcludeTopologies],
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Paint -> OptionValue[FeynArts`Paint]
  ];

  (*Compute the amplitude *)
  ampFA = FeynArts`CreateFeynAmp[diagrams,
    FeynArts`GaugeRules -> OptionValue[FeynArts`GaugeRules],
    FeynArts`PreFactor -> OptionValue[FeynArts`PreFactor],
    FeynArts`Truncated -> OptionValue[FeynArts`Truncated]
  ];

  ampFC = FeynCalc`FCFAConvert[
    ampFA,
    FeynCalc`IncomingMomenta -> OptionValue[FeynCalc`IncomingMomenta],
    FeynCalc`OutgoingMomenta -> OptionValue[FeynCalc`OutgoingMomenta],
    FeynCalc`LorentzIndexNames -> OptionValue[FeynCalc`LorentzIndexNames],
    FeynCalc`LoopMomenta -> OptionValue[FeynCalc`LoopMomenta],
    FeynCalc`UndoChiralSplittings -> OptionValue[FeynCalc`UndoChiralSplittings],
    FeynCalc`SMP -> OptionValue[FeynCalc`SMP],
    FeynCalc`ChangeDimension -> OptionValue[FeynCalc`ChangeDimension],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`TransversePolarizationVectors -> OptionValue[FeynCalc`TransversePolarizationVectors]
  ];

  ampFC = Flatten[ampFC];
  ampFC = If[OptionValue[List], ampFC, Total[ampFC]];
  ampFC = ReplaceAll[ampFC, {FeynArts`FALeviCivita -> FeynCalc`Eps, Global`FAFourVector -> FeynCalc`Momentum}];
  (* subsitute unstable propagator *)
  ampFC = Map[HazmaSubstituteUnstableMediatorPropagator, ampFC];
  ampFC = Map[FeynCalc`PropagatorDenominatorExplicit, ampFC];
  ampFC = Map[FeynCalc`Contract, ampFC];
  ampFC = Map[FeynCalc`DiracSimplify[#, FeynCalc`DiracSubstitute67 -> True]&, ampFC];
  ampFC = Map[FeynCalc`ExpandScalarProduct, ampFC];
  Map[Simplify, ampFC]
];



Options[HazmaComputeAmplitudeSquared] := {
  FeynArts`LoopNumber -> 0,
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`ExcludeTopologies -> {FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`Paint -> False,
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynArts`PreFactor -> -I (2 \[Pi])^(-4 FeynArts`LoopNumber),
  FeynArts`Truncated -> False,
  FeynCalc`ChangeDimension -> 4,
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynCalc`IncomingMomenta -> {},
  List -> False,
  FeynCalc`LoopMomenta -> {},
  FeynCalc`LorentzIndexNames -> {Global`\[Mu], Global`\[Nu], Global`\[Alpha], Global`\[Beta], Global`\[Rho], Global`\[Sigma]},
  FeynCalc`OutgoingMomenta -> {},
  FeynCalc`SMP -> True,
  FeynCalc`TransversePolarizationVectors -> {},
  FeynCalc`UndoChiralSplittings -> True
};


HazmaComputeAmplitudeSquared[inStates_, outStates_, OptionsPattern[]] := Module[{msqrd, amp, inMomenta, outMomenta},

  (* Create the amplitude *)
  amp = HazmaComputeAmplitude[
    inStates,
    outStates,
    FeynArts`LoopNumber -> OptionValue[FeynArts`LoopNumber],
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynArts`ExcludeTopologies -> OptionValue[FeynArts`ExcludeTopologies],
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`GaugeRules -> OptionValue[FeynArts`GaugeRules],
    FeynArts`PreFactor -> OptionValue[FeynArts`PreFactor],
    FeynArts`Truncated -> OptionValue[FeynArts`Truncated],
    FeynCalc`ChangeDimension -> OptionValue[FeynCalc`ChangeDimension],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> OptionValue[FeynCalc`IncomingMomenta],
    List -> OptionValue[List],
    FeynCalc`LoopMomenta -> OptionValue[FeynCalc`LoopMomenta],
    FeynCalc`LorentzIndexNames -> OptionValue[FeynCalc`LorentzIndexNames],
    FeynCalc`OutgoingMomenta -> OptionValue[FeynCalc`OutgoingMomenta],
    FeynCalc`SMP -> OptionValue[FeynCalc`SMP],
    FeynCalc`TransversePolarizationVectors -> OptionValue[FeynCalc`TransversePolarizationVectors],
    FeynCalc`UndoChiralSplittings -> OptionValue[FeynCalc`UndoChiralSplittings]
  ];

  FeynCalc`ClearScalarProducts[];
  (* determine incoming/outgoing momenta *)
  inMomenta = OptionValue[FeynCalc`IncomingMomenta];
  outMomenta = OptionValue[FeynCalc`OutgoingMomenta];

  (* Fix momenta - If user didn't specify momenta, populate them based on what FeynCalc does. *)
  For[i = Length[inMomenta] + 1, i <= Length[inStates], i++,
    AppendTo[inMomenta, ToExpression["InMom" <> ToString[i]]]
  ];
  For[i = Length[outMomenta] + 1, i <= Length[outStates], i++,
    AppendTo[outMomenta, ToExpression["OutMom" <> ToString[i]]]
  ];
  (* set masses *)
  For[i = 1, i <= Length[inMomenta], i++,
    FeynCalc`SP[inMomenta[[i]]] = FeynArts`TheMass[inStates[[i]]]^2;
  ];
  For[i = 1, i <= Length[outMomenta], i++,
    FeynCalc`SP[outMomenta[[i]]] = FeynArts`TheMass[outStates[[i]]]^2
  ];

  (* Square it *)
  msqrd = amp * FeynCalc`FCRenameDummyIndices[FeynCalc`ComplexConjugate[amp]];

  (* If any of the in states are fermionic, average over spins *)
  For[i = 1, i <= Length[inStates], i++,
    If[isFermion[inStates[[i]]],
      msqrd = FeynCalc`FermionSpinSum[msqrd, FeynCalc`ExtraFactor -> 1 / 2],
      Continue
    ];
  ];

  (* If any of the out states are fermionic, sum over spins *)
  For[i = 1, i <= Length[outStates], i++,
    If[isFermion[outStates[[i]]],
      msqrd = FeynCalc`FermionSpinSum[msqrd],
      Continue
    ];
  ];

  (* actually perform the fermion spin sums *)
  msqrd = ReplaceAll[msqrd, FeynCalc`DiracTrace -> FeynCalc`TR];

  (* If any of the states are vectors, do a sum over polarizations. *)
  (* Add 1/2 (1/3) massless (massive) inital state vectrors *)
  For[i = 1, i <= Length[inStates], i++,
    If[isVector[inStates[[i]]],
      If[FeynArts`TheMass[inStates[[i]]] === 0,
        msqrd = (1 / 2)FeynCalc`DoPolarizationSums[msqrd, inMomenta[[i]], 0],
        msqrd = (1 / 3)FeynCalc`DoPolarizationSums[msqrd, inMomenta[[i]]];
      ],
      Continue
    ];
  ];

  For[i = 1, i <= Length[outStates], i++,
    If[isVector[outStates[[i]]],
      If[FeynArts`TheMass[outStates[[i]]] === 0,
        msqrd = FeynCalc`DoPolarizationSums[msqrd, outMomenta[[i]], 0],
        msqrd = FeynCalc`DoPolarizationSums[msqrd, outMomenta[[i]]];
      ],
      Continue
    ];
  ];

  msqrd = FeynCalc`Contract[msqrd];
  msqrd = Simplify[msqrd];
  FeynCalc`ClearScalarProducts[];
  msqrd
];


End[]

EndPackage[]