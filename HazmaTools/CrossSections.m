(* Mathematica Package *)

(* :Title: Cross Sections *)
(* :Context: HazmaTools` *)
(* :Author: Adam Coogan and Logan Morrison *)
(* :Date: 2019-11-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 Adam Coogan and Logan Morrison *)
(* :Keywords: Hazma, Quantum Field Theory, FeynCalc, FeynArts *)
(* :Discussion: File for computing cross-sections in the Hazma framework *)

BeginPackage["HazmaTools`"];

HazmaComputeCrossSection22::usage = "HazmaComputeCrossSection22[]";
HazmaComputeCrossSection22TimesVelocity::usage = "HazmaComputeCrossSection22TimesVelocity[]";
(* Needs work.
  HazmaComputeCrossSection::usage = "HazmaComputeCrossSection[]";
*)

Begin["`Private`"];

IntegralTable[t_, tmin_, tmax_] := {(*** Integral Table ***)
  (** Integrals of t^n / X **)
  (* t^n / (t+a)(t+b) and n >= 1 *)
  c_. t^n_ / ((x_.t + a_) * (y_.t + b_)) :> c / (x * y) * ((
    ((-1)^n * (a / x)^n) / (b / y - a / x)) * Log[(tmax + a / x) / (tmin + a / x)] -
      (((-1)^n * (b / y)^n) / (b / y - a / x)) *
          Log[(tmax + b / y) / (tmin + b / y)] + Sum[Binomial[n, k] * (-1)^k *
      (((a / x)^k * (tmax + a / x)^(n - k) - (b / y)^k * (tmax + b / y)^(n - k)) / ((b / y - a / x) * (n - k)) -
          ((a / x)^k * (tmin + a / x)^(n - k) - (b / y)^k * (tmin + b / y)^(n - k)) / ((b / y - a / x) * (n - k))),
    {k, 0, n - 1}]) /; n > 0 && a / x =!= b / y && FreeQ[c, t] && UnsameQ[x, 0] && UnsameQ[y, 0],
  (* t^n / (t+a)^2 and n >= 2 *)
  c_. t^n_ / (x_.t + a_)^2 :> c / x^2 * ((-1)^(n + 1) * (a / x)^n * (1 / (tmax + a / x) -
      1 / (tmin + a / x)) + n * (-1)^(n - 1) * (a / x)^(n - 1) *
      Log[(tmax + a / x) / (tmin + a / x)] + Sum[Binomial[n, k] * (-1)^k * (a / x)^k *
      (((tmax + a / x)^(n - k - 1) - (tmin + a / x)^(n - k - 1)) / (n - k - 1)),
    {k, 0, n - 2}]) /; n > 1 && FreeQ[c, t] && UnsameQ[x, 0],
  (* t^n / (t+a) and n >= 1 *)
  c_.t^n_ / (x_.t + a_) :> c / x * ((-(a / x))^n * Log[(tmax + a / x) / (tmin + a / x)] +
      Sum[Binomial[n, k] * (-(a / x))^k * (((tmax + a / x)^(n - k) - (tmin + a / x)^(n - k)) / (n - k)),
        {k, 0, n - 1}]) /; n > 0 && FreeQ[c, t] && UnsameQ[x, 0],

  (** Special Case: n = 1 **)
  (* t / (t+a)(t+b) *)
  c_. t / ((x_.t + a_) * (y_.t + b_)) :> c / (x * y) * ((-(a / x) / (b / y - a / x)) * Log[(tmax + a / x) / (tmin + a / x)] -
      (-(b / y) / (b / y - a / x)) * Log[(tmax + b / y) / (tmin + b / y)]) /;
      a / x =!= b / y && FreeQ[c, t] && UnsameQ[x, 0] && UnsameQ[y, 0],
  (* t / (t+a)^2 **)
  c_. t / (x_.t + a_)^2 :> c / x * ((a / x) * (1 / (tmax + a / x) -
      1 / (tmin + a / x)) + Log[(tmax + a / x) / (tmin + a / x)]) /; FreeQ[c, t] && UnsameQ[x, 0],
  (* t / (t+a) *)
  c_. t / (x_.t + a_) :> c / x * (-(a / x)) * Log[(tmax + a / x) / (tmin + a / x)] /; FreeQ[c, t] && UnsameQ[x, 0],

  (** Special Case: n = 0 **)
  (* 1 / (t+a)(t+b) *)
  c_. 1 / ((x_.t + a_) * (y_.t + b_)) :> c / (x * y) * ((1 / (b / y - a / x)) * Log[(tmax + a / x) / (tmin + a / x)] -
      (1 / (b / y - a / x)) * Log[(tmax + b / y) / (tmin + b / y)]) /; FreeQ[c, t] && UnsameQ[x, 0] && UnsameQ[y, 0],
  (* 1 / (t+a)^2 *)
  c_. 1 / (x_.t + a_)^2 :> c / x^2 * (-(1 / (tmax + a / x) - 1 / (tmin + a / x))) /; FreeQ[c, t] && UnsameQ[x, 0],
  (* 1 / (t+a) *)
  c_. 1 / (x_.t + a_) :> c / x * Log[(tmax + a / x) / (tmin + a / x)] /; FreeQ[c, t] && UnsameQ[x, 0],

  (* No match, just integrate *)
  a__ :> Integrate[a, {t, tmin, tmax}, GenerateConditions -> False]
};

CrossSectionIntegrate[exp_, {t_, tmin_, tmax_}] := Module[{integral},
  (* expand out and group terms into list *)
  integral = List @@ Expand[exp];
  (* factor the denominators *)
  integral = Map[Factor[#, GaussianIntegers -> True] &, integral];
  (* replace terms we know with their integrations *)
  integral = Map[ReplaceAll[#, IntegralTable[t, tmin, tmax]] &, integral];
  Together[Total[integral]]
];

Options[HazmaComputeCrossSection22] := {
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynArts`Paint -> False,
  FeynArts`ColumnsXRows -> OptionValue[FeynArts`Paint, FeynArts`ColumnsXRows],
  FeynCalc`FinalSubstitutions -> {},
  FeynCalc`SMP -> True
};

HazmaComputeCrossSection22[inStates_, outStates_, Q_, OptionsPattern[]] := Module[{msqrd, p1, p2, p3, p4, m1, m2, m3, m4, s, t, u, E1, E3, pi, pf, preFactor, tmin, tmax},
  ClearScalarProducts[];
  (* Create the amplitude *)
  msqrd = HazmaComputeAmplitudeSquared[
    inStates,
    outStates,
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows],
    FeynArts`GaugeRules -> OptionValue[FeynArts`GaugeRules],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {p1, p2},
    FeynCalc`OutgoingMomenta -> {p3, p4},
    FeynCalc`SMP -> OptionValue[FeynCalc`SMP]
  ];

  m1 = FeynArts`TheMass[inStates[[1]]];
  m2 = FeynArts`TheMass[inStates[[2]]];
  m3 = FeynArts`TheMass[outStates[[1]]];
  m4 = FeynArts`TheMass[outStates[[2]]];

  FeynCalc`SetMandelstam[s, t, u, p1, p2, -p3, -p4, m1, m2, m3, m4];
  msqrd = Simplify[msqrd];
  msqrd = ReplaceAll[msqrd, {u -> m1^2 + m2^2 + m3^2 + m4^2 - s - t}];

  E1 = (s + m1^2 - m2^2) / (2 * Sqrt[s]);
  E3 = (s + m3^2 - m4^2) / (2 * Sqrt[s]);
  pi = Sqrt[E1^2 - m1^2]; (* mag. of initial state 3-momentum *)
  pf = Sqrt[E3^2 - m3^2]; (* mag. of final state 3-momentum *)
  (*msqrd=Simplify[ReplaceAll[msqrd,{t->m1^2+m3^2-2(E1*E3-z p1*p3)}]];*)

  preFactor = Simplify[(2 * \[Pi]) / (64 * \[Pi]^2 * Q^2) * (pf / pi) * (1 / (2 * pi * pf))];
  (* If final state particles are identical, divide by 2 *)
  If[outStates[[1]] === outStates[[2]], preFactor = preFactor / 2, Continue];

  tmax = Simplify[m1^2 + m3^2 - 2E1 * E3 + 2 * pi * pf];
  tmin = Simplify[m1^2 + m3^2 - 2E1 * E3 - 2 * pi * pf];


  ReplaceAll[preFactor * Integrate[msqrd, {t, tmin, tmax}, GenerateConditions -> False], {s -> Q^2}]
];

Options[HazmaComputeCrossSection22TimesVelocity] := {
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`ExcludeTopologies -> {FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynArts`Paint -> False,
  FeynArts`ColumnsXRows -> OptionValue[FeynArts`Paint, FeynArts`ColumnsXRows],
  FeynCalc`ChangeDimension -> 4,
  FeynCalc`FinalSubstitutions -> {},
  FeynCalc`SMP -> True
};

HazmaComputeCrossSection22TimesVelocity[inState_, outStates_, OptionsPattern[]] := Module[{msqrd, p1, p2, p3, p4, m, m3, m4, s, t, u, E1, E3, preFactor, z, \[Beta]f},
  If[Length[inState] === 1, Continue, Throw[$Failed, HazmaComputeCrossSection22TimesVelocity::InvalidInState]];
  ClearScalarProducts[];
  (* Create the amplitude *)
  msqrd = HazmaComputeAmplitudeSquared[
    {inState, -inState},
    outStates,
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynArts`ExcludeTopologies -> OptionValue[FeynArts`ExcludeTopologies],
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows],
    FeynArts`GaugeRules -> OptionValue[FeynArts`GaugeRules],
    FeynCalc`ChangeDimension -> OptionValue[FeynCalc`ChangeDimension],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {p1, p2},
    FeynCalc`OutgoingMomenta -> {p3, p4},
    FeynCalc`SMP -> OptionValue[FeynCalc`SMP]
  ];

  m = FeynArts`TheMass[inState];
  m3 = FeynArts`TheMass[outStates[[1]]];
  m4 = FeynArts`TheMass[outStates[[2]]];

  FeynCalc`SetMandelstam[s, t, u, p1, p2, -p3, -p4, m, m, m3, m4];
  msqrd = Simplify[msqrd];
  msqrd = ReplaceAll[msqrd, {u -> m^2 + m^2 + m3^2 + m4^2 - s - t}];

  (* Assume #1 is at rest: p1 = (m,0,0,0)*)
  (* Then, p2 = (E2,0,0,p) where E2 = s/(2m)-m *)
  E1 = m;
  E3 = (m3^2 - m4^2 + s - Sqrt[-4 * m^2 + s] * Sqrt[-2 * (m3^2 + m4^2) + (m3^2 - m4^2)^2 / s + s] * z) / (4 * m);
  (*msqrd=Simplify[ReplaceAll[msqrd,{t->m1^2+m3^2-2(E1*E3-z p1*p3)}]];*)

  (* Gondolo-gelmini Eqn.(3.24) *)
  \[Beta]f = Sqrt[(1 - (m3 + m4)^2 / s)(1 - (m3 - m4)^2 / s)];
  preFactor = Simplify[(2 * \[Pi]) / (64 * \[Pi]^2 * (s - 2 * m^2)) * \[Beta]f];
  (* If final state particles are identical, divide by 2 *)
  If[outStates[[1]] === outStates[[2]], preFactor = preFactor / 2, Continue];

  msqrd = ReplaceAll[msqrd, {t -> m^2 + (Sqrt[-4 * m^2 + s] * (m3^2 - m4^2 + s) *
      (-(s * Sqrt[(m3^4 + (m4^2 - s)^2 - 2 * m3^2 * (m4^2 + s)) / s]) +
          Sqrt[s * (m3^4 + (m4^2 - s)^2 - 2 * m3^2 * (m4^2 + s))]) * z) / (8 * m^2 * s) +
      (m3^2 + m4^2 - s + Sqrt[-4 * m^2 + s] * Sqrt[-2 * (m3^2 + m4^2) +
          (m3^2 - m4^2)^2 / s + s] * z) / 2}];

  ReplaceAll[preFactor * Integrate[msqrd, {z, -1, 1}, GenerateConditions -> False], {s -> 4 * m^2 (1 + Global`\[Epsilon])}]

];

Options[HazmaComputeCrossSection] := {
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`GaugeRules -> _FeynArts`FAGaugeXi -> 1,
  FeynArts`Paint -> False,
  FeynArts`ColumnsXRows -> OptionValue[FeynArts`Paint, FeynArts`ColumnsXRows],
  FeynCalc`FinalSubstitutions -> {},
  FeynCalc`SMP -> True
};

HazmaComputeCrossSection[inStates_, outStates_, Q_, OptionsPattern[]] := Module[{msqrd, p1, p2, p3, p4, m1, m2, m3, m4, s, t, u, preFactor, tmin, tmax},
  ClearScalarProducts[];
  (* Create the amplitude *)
  msqrd = HazmaComputeAmplitudeSquared[
    inStates,
    outStates,
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles],
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Paint -> OptionValue[FeynArts`Paint],
    FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows],
    FeynArts`GaugeRules -> OptionValue[FeynArts`GaugeRules],
    FeynCalc`FinalSubstitutions -> OptionValue[FeynCalc`FinalSubstitutions],
    FeynCalc`IncomingMomenta -> {p1, p2},
    FeynCalc`OutgoingMomenta -> {p3, p4},
    FeynCalc`SMP -> OptionValue[FeynCalc`SMP]
  ];

  m1 = FeynArts`TheMass[inStates[[1]]];
  m2 = FeynArts`TheMass[inStates[[2]]];
  m3 = FeynArts`TheMass[outStates[[1]]];
  m4 = FeynArts`TheMass[outStates[[2]]];

  FeynCalc`SetMandelstam[s, t, u, p1, p2, -p3, -p4, m1, m2, m3, m4];
  msqrd = Simplify[msqrd];
  msqrd = ReplaceAll[msqrd, {u -> m1^2 + m2^2 + m3^2 + m4^2 - s - t}];

  preFactor = 1 / (16 * Pi * KallenLambda[s, m1^2, m2^2]);
  (* If final state particles are identical, divide by 2 *)
  If[outStates[[1]] === outStates[[2]], preFactor = preFactor / 2, Continue];

  tmin = m1^2 + m3^2 - (Sqrt[KallenLambda[s, m3^2, m4^2] * KallenLambda[s, m1^2, m2^2]] + (s + m1^2 - m2^2) * (s + m3^2 - m4^2)) / (2 * s);
  tmax = m1^2 + m3^2 + (Sqrt[KallenLambda[s, m3^2, m4^2] * KallenLambda[s, m1^2, m2^2]] - (s + m1^2 - m2^2) * (s + m3^2 - m4^2)) / (2 * s);

  ReplaceAll[preFactor * CrossSectionIntegrate[msqrd, {t, tmin, tmax}], {s -> Q^2}]
];

End[];

EndPackage[]
