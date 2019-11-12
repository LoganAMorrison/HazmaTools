(* Mathematica Package *)

(* :Title: Widths *)
(* :Context: HazmaTools` *)
(* :Author: Adam Coogan and Logan Morrison *)
(* :Date: 2019-11-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 Adam Coogan and Logan Morrison *)
(* :Keywords: Hazma, Quantum Field Theory, FeynCalc, FeynArts *)
(* :Discussion: File for computing decay widths in the Hazma framework *)

BeginPackage["HazmaTools`"]

HazmaComputeWidth::usage = "HazmaComputeWidth[inState, outStates] Compute the partial decay width \
of inState -> outState. Note inState must be of length 1 and outState must be of length 2 or 3."


Begin["`Private`"]

HazmaComputeWidth::InvalidInStates = "Invalid number of inStates passed to HazmaComputeWidth. Use \
1 state."

HazmaComputeWidth::InvalidOutStates = "Invalid number of outStates passed to HazmaComputeWidth. Use \
2 or 3 states."

Options[HazmaComputeWidth] = {
  FeynArts`LoopNumber -> 0,
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`ExcludeTopologies -> {FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`Paint -> False,
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynCalc`ChangeDimension -> 4,
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynCalc`LoopMomenta -> {},
  FeynCalc`SMP -> True
};

HazmaComputeWidth[inState_, outStates_, opt : OptionsPattern[]] := Module[{},
  If[Length[inState] != 1,
    Message[HazmaComputeWidth::InvalidInStates];
    Throw[$Failed]
  ];
  If[Length[outStates] == 2, Return[HazmaComputeWidth12[inState, outStates, opt]]];
  If[Length[outStates] == 3, Return[HazmaComputeWidth13[inState, outStates, opt]]];
  Message[HazmaComputeWidth::InvalidOutStates];
  Throw[$Failed]
];




Options[HazmaComputeWidth12] = {
  FeynArts`LoopNumber -> 0,
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`ExcludeTopologies -> {FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`Paint -> False,
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynCalc`ChangeDimension -> 4,
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynCalc`LoopMomenta -> {},
  FeynCalc`SMP -> True
};

HazmaComputeWidth12[inState_, outStates_, OptionsPattern[]] := Module[{msqrd, P, p1, p2, M, m1, m2, pf, preFactor, E1, E2},
  FeynCalc`ClearScalarProducts[];
  (* Create the amplitude *)
  msqrd = HazmaComputeAmplitudeSquared[
    {inState},
    outStates,
    FeynArts`LoopNumber -> OptionValue[FeynArts`LoopNumber],
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynArts`ExcludeTopologies -> OptionValue[FeynArts`ExcludeTopologies],
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`GaugeRules -> OptionValue[FeynArts`GaugeRules],
    FeynCalc`ChangeDimension -> OptionValue[FeynCalc`ChangeDimension],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {P},
    FeynCalc`LoopMomenta -> OptionValue[FeynCalc`LoopMomenta],
    FeynCalc`OutgoingMomenta -> {p1, p2},
    FeynCalc`SMP -> OptionValue[FeynCalc`SMP]
  ];

  M = FeynArts`TheMass[inState];
  m1 = FeynArts`TheMass[outStates[[1]]];
  m2 = FeynArts`TheMass[outStates[[2]]];

  E1 = (M^2 + m1^2 - m2^2) / (2 * M);
  E2 = (M^2 - m1^2 + m2^2) / (2 * M);

  FeynCalc`SP[P, P] = M^2;
  FeynCalc`SP[p1, p1] = m1^2;
  FeynCalc`SP[p2, p2] = m2^2;
  FeynCalc`SP[P, p1] = M * E1;
  FeynCalc`SP[P, p2] = M * E2;
  FeynCalc`SP[p1, p2] = (M^2 - m1^2 - m2^2) / 2;

  msqrd = Simplify[msqrd];

  pf = Sqrt[E1^2 - m1^2]; (* mag. of final state 3-momentum *)

  preFactor = Simplify[(4 * \[Pi]) / (2 * M) * pf / (16 * Pi^2 * M)];
  (* If final state particles are identical, divide by 2 *)
  preFactor *= If[outStates[[1]] === outStates[[2]], 1 / 2, 1];
  preFactor * msqrd

];


Options[HazmaComputeWidth13] := {
  FeynArts`LoopNumber -> 0,
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`ExcludeTopologies -> {FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`Paint -> False,
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynCalc`ChangeDimension -> 4,
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynCalc`LoopMomenta -> {},
  FeynCalc`SMP -> True
};

HazmaComputeWidth13[inState_, outStates_, OptionsPattern[]] := Module[{msqrd, P, p1, p2, p3, M, m1, m2, m3, preFactor, tbounds, tintegral, sintegral},
  ClearScalarProducts[];
  (* Create the amplitude *)
  msqrd = HazmaComputeAmplitudeSquared[
    {inState},
    outStates,
    FeynArts`LoopNumber -> OptionValue[FeynArts`LoopNumber],
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynArts`ExcludeTopologies -> OptionValue[FeynArts`ExcludeTopologies],
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`GaugeRules -> OptionValue[FeynArts`GaugeRules],
    FeynCalc`ChangeDimension -> OptionValue[FeynCalc`ChangeDimension],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {P},
    FeynCalc`LoopMomenta -> OptionValue[FeynCalc`LoopMomenta],
    FeynCalc`OutgoingMomenta -> {p1, p2, p3},
    FeynCalc`SMP -> OptionValue[FeynCalc`SMP]
  ];

  M = FeynArts`TheMass[inState];
  m1 = FeynArts`TheMass[outStates[[1]]];
  m2 = FeynArts`TheMass[outStates[[2]]];
  m3 = FeynArts`TheMass[outStates[[3]]];

  FeynCalc`SetMandelstam[s, t, u, P, -p1, -p2, -p3, M, m1, m2, m3];
  msqrd = Simplify[msqrd];
  msqrd = ReplaceAll[msqrd, u -> {M^2 + m1^2 + m2^2 + m3^2 - s - t}];

  tbounds = ComputeMandelstamTBound[M, m1, m2, m3, s, t];

  (* Integrate over t *)
  tintegral = Integrate[msqrd, t, GenerateConditions -> False];
  tintegral = ReplaceAll[tintegral, tbounds[[2]]] - ReplaceAll[tintegral, tbounds[[1]]];
  (* Integrate over s *)
  sintegral = Integrate[tintegral, s, GenerateConditions -> False];
  sintegral = ReplaceAll[sintegral, s -> (M - m1)^2] - ReplaceAll[sintegral, s -> (m2 + m3)^2];

  preFactor = Simplify[1 / (2 * M) * 1 / (16 * M^2 * (2 * Pi)^3)];
  (* If final state particles are identical, divide by 2 *)
  preFactor /= Factorial[4 - CountDistinct[outStates]];
  preFactor * sintegral

];

End[]


EndPackage[]