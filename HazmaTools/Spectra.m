(* Mathematica Package *)

(* :Title: Spectra *)
(* :Context: HazmaTools` *)
(* :Author: Adam Coogan and Logan Morrison *)
(* :Date: 2019-11-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 Adam Coogan and Logan Morrison *)
(* :Keywords: Hazma, Quantum Field Theory, FeynCalc, FeynArts *)
(* :Discussion: File for computing gamma-ray spectra in the Hazma framework *)

BeginPackage["HazmaTools`"];

HazmaComputeDNDE::usage = "HazmaComputeDNDE[inStates, outStates, Q] or HazmaComputeDNDE[outStates] Compute the \
differential radiative gamma-ray spectrum from inStates->outStates or mediator->outStates. Note that outStates should \
not include Photon.";


Begin["`Private`"];

HazmaComputeDNDE::InvalidModel = "Implementation for model `1 not yet implemented";

HazmaComputeDNDE::InvaliedOutStates = "Invalid outStates. Must be of length 2. \
Do not include photon in list.";

Options[HazmaComputeDNDE] := {
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`Paint -> False,
  FeynArts`ColumnsXRows -> OptionValue[FeynArts`Paint, FeynArts`ColumnsXRows],
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynArts`ExcludeParticles -> {Photon}
};

HazmaComputeDNDE[inStates_, outStates_, Q_, opt : OptionsPattern[]] := Module[{},
  If[$HazmaModel === "scalar", Return[ScalarMediatorComputeDNDE[inStates, outStates, Q, opt]]];
  If[$HazmaModel === "vector", Return[VectorMediatorComputeDNDE[inStates, outStates, Q, opt]]];
  Message[HazmaComputeDNDE::InvalidModel, $HazmaModel];
  Throw[$Failed]
];


HazmaComputeDNDE[outStates_, OptionsPattern[]] := Module[{mediator, msqrd, width, P, p1, p2, p3, M, m1, m2, m3, s, t, u, tub, tlb, tintegral, dndE, preFactor},
  If[Length[outStates] != 2,
    Message[HazmaComputeDNDE::InvaliedOutStates];
    Throw[$Failed]
  ];
  If[$HazmaModel === "scalar", mediator = ScalarMediator];
  If[$HazmaModel === "vector", mediator = VectorMediator];

  (* Compute tree-level width *)
  width = HazmaComputeWidth[
    mediator,
    outStates,
    FeynArts`Adjacencies -> OptionValue[Adjacencies],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles]
  ];

  (* Determine masses *)
  M = FeynArts`TheMass[mediator];
  m1 = 0;
  m2 = FeynArts`TheMass[outStates[[1]]];
  m3 = FeynArts`TheMass[outStates[[2]]];

  msqrd = HazmaComputeAmplitudeSquared[
    {mediator},
    {Photon, outStates[[1]], outStates[[2]]},
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows],
    (* ignore diagrams with internal photon. They are higher order in \[Alpha] *)
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynCalc`ChangeDimension -> 4,
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {P},
    FeynCalc`OutgoingMomenta -> {p1, p2, p3},
    FeynCalc`SMP -> True
  ];
  FeynCalc`SetMandelstam[s, t, u, P, -p1, -p2, -p3, M, m1, m2, m3];

  msqrd = ReplaceAll[msqrd, u -> M^2 + m1^2 + m2^2 + m3^2 - s - t];

  {tub, tlb} = ComputeMandelstamTBound[M, m1, m2, m3, s, t];
  tintegral = Integrate[msqrd, t, GenerateConditions -> False];
  tintegral = Simplify[(tintegral /. tub) - (tintegral /. tlb), Assumptions -> {M > 0, M > m2 + m3 > 0}];

  preFactor = -1 / (16 * M^2 * (2 * Pi)^3) * 1 / width;
  dndE = tintegral * preFactor;

  dndE = ReplaceAll[dndE, {s -> M^2 - 2 * M * Global`E\[Gamma], Global`qe -> Sqrt[4 * Pi * Global`alphaEM]}];
  Simplify[dndE, Assumptions -> {M > 0, M - (m2 + m3)^2 / M > 2 * Global`E\[Gamma] > 0, m1 >= 0, m2 >= 0, m3 >= 0}]
];


Options[ScalarMediatorComputeDNDE] := {
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`Paint -> False,
  FeynArts`ColumnsXRows -> OptionValue[FeynArts`Paint, FeynArts`ColumnsXRows],
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynArts`ExcludeParticles -> {Photon}
};

ScalarMediatorComputeDNDE[inStates_, outStates_, Q_, OptionsPattern[]] := Module[{msqrd, \[Sigma]0, p1, p2, p3, p4, k, m1, m2, m3, m4, s, t, u, p, tub, tlb, tintegral, dndE, preFactor},
  If[Length[inStates] == 2, Continue, Throw[$Failed, ScalarMediatorComputeDNDE::InvalidInStates]];
  If[Length[outStates] == 2, Continue, Throw[$Failed, ScalarMediatorComputeDNDE::InvalidOutStates]];
  (* Compute tree-level cross section *)
  \[Sigma]0 = HazmaComputeCrossSection22[
    inStates,
    outStates,
    Q,
    Adjacencies -> OptionValue[Adjacencies],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles]
  ];
  (* Determine masses *)
  m1 = FeynArts`TheMass[inStates[[1]]];
  m2 = FeynArts`TheMass[inStates[[2]]];
  m3 = FeynArts`TheMass[outStates[[1]]];
  m4 = FeynArts`TheMass[outStates[[2]]];
  (* Set masses of four-momenta *)
  FeynCalc`SP[p1, p1] = m1^2;
  FeynCalc`SP[p2, p2] = m2^2;
  FeynCalc`SP[p3, p3] = m3^2;
  FeynCalc`SP[p4, p4] = m4^2;
  FeynCalc`SP[k, k] = 0;
  msqrd = HazmaComputeAmplitudeSquared[
    inStates,
    {outStates[[1]], outStates[[2]], Photon},
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows],
    (* ignore diagrams with internal photon. They are higher order in \[Alpha] *)
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynCalc`ChangeDimension -> 4,
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {p1, p2},
    FeynCalc`OutgoingMomenta -> {p3, p4, k},
    FeynCalc`SMP -> True
  ];
  (*
  Set the Mandelstam variables:
    s = (P - k)^2 = (p3 + p4)^2 = m3^2 + m4^2 + 2p3.p4
      => p3.p4 = (1/2)(s - m3^2 - m4^2)
    t = (P - p3)^2 = (p4 + k)^2 = m4^2 + 2k.p4
      => k.p4 = (1/2)(t - m4^2)
    u = (P - p4)^2 = (p3 + k)^2 = m3^2 + 2p3.k
      => k.p3 = (1/2)(u - m3^2)
    Q^2 = (p1 + p2)^2 = m1^2 + m2^2 + 2p1.p2
      => p1.p2 = (1/2)(Q^2 - m1^2 - m2^2)
  *)
  FeynCalc`SP[p3, p4] = (1 / 2) * (s - m3^2 - m4^2);
  FeynCalc`SP[k, p4] = (1 / 2) * (t - m4^2);
  FeynCalc`SP[k, p3] = (1 / 2) * (u - m3^2);
  FeynCalc`SP[p1, p2] = (1 / 2) * (Q^2 - m1^2 - m2^2);

  (* Replace u: s + t + u = Q^2 + m3^2 + m4^2 + 0 *)
  msqrd = ReplaceAll[msqrd, u -> Q^2 + m3^2 + m4^2 - s - t];

  (* Integrate over t *)
  {tub, tlb} = ComputeMandelstamTBound[Q, 0, m3, m4, s, t];
  tintegral = Integrate[msqrd, t, GenerateConditions -> False];
  tintegral = Simplify[(tintegral /. tub) - (tintegral /. tlb), Assumptions -> {Q > 0, Q > m1 + m2 > 0, Q > m3 + m4 > 0}];

  (*
  Compute the prefactor of dNdE:
  The differential cross section we have is in terms of s and t:
    \[DifferentialD]\[Sigma] = 1/(2E1*2E2*vrel) * M^2 * \[DifferentialD]\[CapitalPi]
  where \[DifferentialD]\[CapitalPi] is:
    \[DifferentialD]\[CapitalPi] = \[DifferentialD]s*\[DifferentialD]t / (16*Q^2(2*Pi)^3)
  Using:
    vrel = p/E1 + p/E2 = p(E1+E2)/(E1*E2) = p Q/(E1*E2)
  and
    p = Sqrt[E1^2-m1^2], E1=(Q^2+m1^2-m2^2)/(2Q);
  we find that:
    \[DifferentialD]\[Sigma] = 1/(4*p*Q) * M^2 * \[DifferentialD]s*\[DifferentialD]t / (16*Q^2(2*Pi)^3)
  Integrating over t, we find:
    \[DifferentialD]\[Sigma]/\[DifferentialD]s = 1/(4*p*Q*16*Q^2(2*Pi)^3)\[Integral] \[DifferentialD]t*M^2
  To convert to \[DifferentialD]\[Sigma]/\[DifferentialD]E, we use:
    \[DifferentialD]s/\[DifferentialD]E = \[DifferentialD]/\[DifferentialD]E (Q^2-2QE) = -2Q
  Thus:
    \[DifferentialD]\[Sigma]/\[DifferentialD]E = (-2Q)/(4*p*Q*16*Q^2(2*Pi)^3)\[Integral] \[DifferentialD]t*M^2
  *)
  p = Sqrt[(m1 - m2 - Q) * (m1 + m2 - Q) * (m1 - m2 + Q) * (m1 + m2 + Q)] / (2 * Q);
  preFactor = -(2 * Q) / (4 * Q * p) * 1 / (16 * Q^2 * (2 * Pi)^3) * 1 / \[Sigma]0;
  dndE = tintegral * preFactor;

  (* Convert s to E => s = (Q^2-2*Q*E) *)
  dndE = ReplaceAll[dndE, {s -> Q^2 - 2 * Q * Global`E\[Gamma], Global`qe -> Sqrt[4 * Pi * Global`alphaEM]}];
  Simplify[dndE, Assumptions -> {Q > 0, Q - (m3 + m4)^2 / Q > 2 * Global`E\[Gamma] > 0, m1 > 0, m2 > 0, m3 > 0, m4 > 0}]
];




Options[VectorMediatorComputeDNDE] := {
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`Paint -> False,
  FeynArts`ColumnsXRows -> OptionValue[FeynArts`Paint, FeynArts`ColumnsXRows],
  FeynCalc`FinalSubstitutions -> If[MemberQ[{"scalar", "vector"}, $HazmaModel], {Global`M$FACouplings}, {}],
  FeynArts`ExcludeParticles -> {Photon}
};

VectorMediatorComputeDNDE[inStates_, outStates_, Q_, OptionsPattern[]] := Module[{newOutStates, msqrd, \[Sigma]0, P, p1, p2, p3, p4, k, m1, m2, m3, m4, s, t, u, p, tub, tlb, tintegral, dndE, preFactor, X\[Mu]\[Nu], L\[Mu]\[Nu], ampFS, ampIS, X, L},
  \[Sigma]0 = HazmaComputeCrossSection22[inStates, outStates, Q, Adjacencies -> OptionValue[Adjacencies]];
  newOutStates = outStates;
  AppendTo[newOutStates, state\[Gamma]];
  (* Factor the amplitude into L\[Mu]\[Nu]*X\[Mu]\[Nu]/propagator *)
  (* Initial state M(xx\[Rule]V)^2/4 *)
  ampIS = FeynCalc`Contract[HazmaComputeAmplitude[
    inStates,
    {statev},
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynCalc`ChangeDimension -> 4,
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {p1, p2},
    List -> False,
    FeynCalc`OutgoingMomenta -> {P},
    FeynCalc`SMP -> True
  ]];
  (* Final state piece *)
  ampFS = FeynCalc`Contract[HazmaComputeAmplitude[
    {statev},
    {outStates[[1]], outStates[[2]], Photon},
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynCalc`ChangeDimension -> 4,
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {P},
    List -> False,
    FeynCalc`OutgoingMomenta -> {p3, p4, k},
    FeynCalc`SMP -> True
  ]];
  m1 = FeynArts`TheMass[inStates[[1]]];
  m2 = FeynArts`TheMass[inStates[[2]]];
  m3 = FeynArts`TheMass[outStates[[1]]];
  m4 = FeynArts`TheMass[outStates[[2]]];
  (* set masses *)
  FeynCalc`SP[p1, p1] = m1^2;
  FeynCalc`SP[p2, p2] = m2^2;
  FeynCalc`SP[p3, p3] = m3^2;
  FeynCalc`SP[p4, p4] = m4^2;
  FeynCalc`SP[k, k] = 0;

  (* strip off polarization vectors *)
  ampIS = ReplaceAll[ampIS, {
    FeynCalc`DiracGamma[FeynCalc`Momentum[FeynCalc`Polarization[P, -I]]] :> FeynCalc`GA[\[Mu]],
    FeynCalc`Pair[a__FeynCalc`Momentum, FeynCalc`Momentum[FeynCalc`Polarization[P, I]]] :> FeynCalc`Pair[a, FeynCalc`LorentzIndex[\[Mu]]],
    FeynCalc`Pair[FeynCalc`Momentum[FeynCalc`Polarization[P, I]], a__FeynCalc`Momentum] :> FeynCalc`Pair[a, FeynCalc`LorentzIndex[\[Mu]]]}
  ];
  (* square amplitude and replace mu with new Lorentz index*)
  L\[Mu]\[Nu] = ampIS * ReplaceAll[FeynCalc`ComplexConjugate[ampIS], \[Mu] -> \[Nu]];
  (* average over initial state fermions spins *)
  L\[Mu]\[Nu] = ReplaceAll[FeynCalc`FermionSpinSum[L\[Mu]\[Nu], FeynCalc`ExtraFactor -> 1 / 4], FeynCalc`DiracTrace -> FeynCalc`TR];
  (* Compute L using L\[Mu]\[Nu] = (P^\[Mu]P^\[Nu]-g^\[Mu]\[Nu]P^2)L => L = -g^\[Mu]\[Nu]Subscript[L, \[Mu]\[Nu]]/3Q^2 *)
  L = Simplify[FeynCalc`Contract[L\[Mu]\[Nu] * FeynCalc`MT[\[Mu], \[Nu]]]] / (3Q^2);

  (* Final state Integrate[M(V\[Rule]finalstate)^2/4,phasespace] *)

  (* strip off polarization vectors *)
  ampFS = ReplaceAll[ampFS, {
    FeynCalc`DiracGamma[FeynCalc`Momentum[FeynCalc`Polarization[P, I]]] :> FeynCalc`GA[\[Mu]],
    FeynCalc`Pair[a__FeynCalc`Momentum, FeynCalc`Momentum[FeynCalc`Polarization[P, I]]] :> FeynCalc`Pair[a, FeynCalc`LorentzIndex[\[Mu]]],
    FeynCalc`Pair[FeynCalc`Momentum[FeynCalc`Polarization[P, I]], a__FeynCalc`Momentum] :> FeynCalc`Pair[a, FeynCalc`LorentzIndex[\[Mu]]]}
  ];
  (* square amplitude and replace mu with new Lorentz index*)
  X\[Mu]\[Nu] = ampFS * ReplaceAll[FeynCalc`ComplexConjugate[ampFS], \[Mu] -> \[Nu]];
  (* sum over final state fermions spins *)
  X\[Mu]\[Nu] = ReplaceAll[FeynCalc`FermionSpinSum[X\[Mu]\[Nu]], FeynCalc`DiracTrace -> FeynCalc`TR];
  (* sum over final state photon polarizations *)
  X\[Mu]\[Nu] = FeynCalc`DoPolarizationSums[X\[Mu]\[Nu], k, 0];
  (* Compute X using X\[Mu]\[Nu] = (P^\[Mu]P^\[Nu]-g^\[Mu]\[Nu]P^2)X => X = -g^\[Mu]\[Nu]Subscript[X, \[Mu]\[Nu]]/3Q^2 *)
  X = Simplify[FeynCalc`Contract[X\[Mu]\[Nu] * FeynCalc`MT[\[Mu], \[Nu]]]] / (3Q^2);

  (* Compute M^2 = L\[Mu]\[Nu]X\[Mu]\[Nu] = 3Q^2 X L and add vector propagator *)
  msqrd = X * L (3 * Q^4) / ((Q^2 - Global`mv^2)^2 + (Global`mv * Global`widthv)^2);

  (*
  Set the Mandelstam variables:
    s = (P - k)^2 = (p3 + p4)^2 = m3^2 + m4^2 + 2p3.p4
      => p3.p4 = (1/2)(s - m3^2 - m4^2)
    t = (P - p3)^2 = (p4 + k)^2 = m4^2 + 2k.p4
      => k.p4 = (1/2)(t - m4^2)
    u = (P - p4)^2 = (p3 + k)^2 = m3^2 + 2p3.k
      => k.p3 = (1/2)(u - m3^2)
    Q^2 = (p1 + p2)^2 = m1^2 + m2^2 + 2p1.p2
      => p1.p2 = (1/2)(Q^2 - m1^2 - m2^2)
  *)
  FeynCalc`SP[p3, p4] = (1 / 2) * (s - m3^2 - m4^2);
  FeynCalc`SP[k, p4] = (1 / 2) * (t - m4^2);
  FeynCalc`SP[k, p3] = (1 / 2) * (u - m3^2);
  FeynCalc`SP[p1, p2] = (1 / 2) * (Q^2 - m1^2 - m2^2);

  (* Replace u: s + t + u = Q^2 + m3^2 + m4^2 + 0 *)
  msqrd = ReplaceAll[msqrd, u -> Q^2 + m3^2 + m4^2 - s - t];

  (* Integrate over t *)
  {tub, tlb} = ComputeMandelstamTBound[Q, 0, m3, m4, s, t];
  tintegral = Integrate[msqrd, t, GenerateConditions -> False];
  tintegral = Simplify[(tintegral /. tub) - (tintegral /. tlb), Assumptions -> {Q > 0, Q > m1 + m2 > 0, Q > m3 + m4 > 0}];

  (*
  Compute the prefactor of dNdE:
  The differential cross section we have is in terms of s and t:
    \[DifferentialD]\[Sigma] = 1/(2E1*2E2*vrel) * M^2 * \[DifferentialD]\[CapitalPi]
  where \[DifferentialD]\[CapitalPi] is:
    \[DifferentialD]\[CapitalPi] = \[DifferentialD]s*\[DifferentialD]t / (16*Q^2(2*Pi)^3)
  Using:
    vrel = p/E1 + p/E2 = p(E1+E2)/(E1*E2) = p Q/(E1*E2)
  and
    p = Sqrt[E1^2-m1^2], E1=(Q^2+m1^2-m2^2)/(2Q);
  we find that:
    \[DifferentialD]\[Sigma] = 1/(4*p*Q) * M^2 * \[DifferentialD]s*\[DifferentialD]t / (16*Q^2(2*Pi)^3)
  Integrating over t, we find:
    \[DifferentialD]\[Sigma]/\[DifferentialD]s = 1/(4*p*Q*16*Q^2(2*Pi)^3)\[Integral] \[DifferentialD]t*M^2
  To convert to \[DifferentialD]\[Sigma]/\[DifferentialD]E, we use:
    \[DifferentialD]s/\[DifferentialD]E = \[DifferentialD]/\[DifferentialD]E (Q^2-2QE) = -2Q
  Thus:
    \[DifferentialD]\[Sigma]/\[DifferentialD]E = (-2Q)/(4*p*Q*16*Q^2(2*Pi)^3)\[Integral] \[DifferentialD]t*M^2
  *)
  p = Sqrt[(m1 - m2 - Q) * (m1 + m2 - Q) * (m1 - m2 + Q) * (m1 + m2 + Q)] / (2 * Q);
  preFactor = -(2 * Q) / (4 * Q * p) * 1 / (16 * Q^2 * (2 * Pi)^3) * 1 / \[Sigma]0;
  dndE = tintegral * preFactor;

  (* Convert s to E => s = (Q^2-2*Q*E) *)
  dndE = ReplaceAll[dndE, {s -> Q^2 - 2 * Q * Global`E\[Gamma], Global`qe -> Sqrt[4 * Pi * Global`alphaEM]}];
  Simplify[dndE, Assumptions -> {Q > 0, Q - (m3 + m4)^2 / Q > 2 * Global`E\[Gamma] > 0, m1 > 0, m2 > 0, m3 > 0, m4 > 0}]
];



End[];

EndPackage[]