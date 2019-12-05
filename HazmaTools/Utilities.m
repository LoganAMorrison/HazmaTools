(* Mathematica Package *)

(* :Title: Utilities *)
(* :Context: HazmaTools` *)
(* :Author: Adam Coogan and Logan Morrison *)
(* :Date: 2019-11-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 Adam Coogan and Logan Morrison *)
(* :Keywords: Hazma, Quantum Field Theory, FeynCalc, FeynArts *)
(* :Discussion: Utilities for HazmaTools *)

BeginPackage["HazmaTools`"];

$HazmaModel::usage = "$HazmaModel denotes the model which should be used. \
Set `$HazmaModel=\"scalar\"` for scalar mediator, `$HazmaModel=\"vector\"`\ 
for vector mediator, etc. By default, it is set to scalar.";

HazmaSubstituteUnstableMediatorPropagator::usage = "HazmaSubstitutionUnstableMediatorPropagator[exp] \
substitute the mediator propagators with unstable propagators.";

ComputeMandelstamTBound::usage = "ComputeMandelstamTBound[M, m1, m2, m3, s, t] Computes the bounds on \
the Mandelstam variable t = (P-p2)^2 where p2 is the four-momentum associated with m2.";

Begin["`Private`"];

$HazmaModel = "HazmaScalar";

SubstituteUnstableScalarPropagator[exp_] := ReplaceAll[exp, {
  FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum + b__FeynCalc`Momentum, Global`ms]] :> 1 / (FeynCalc`SP[a + b, a + b] - Global`ms^2 + I Global`widths Global`ms),
  FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a_FeynCalc`Momentum - b_FeynCalc`Momentum, Global`ms]] :> 1 / (FeynCalc`SP[a - b, a - b] - Global`ms^2 + I Global`widths Global`ms),
  FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum, Global`ms]] :> 1 / (FeynCalc`SP[a, a] - Global`ms^2 + I Global`widths Global`ms)}
];

SubstituteUnstableVectorPropagator[exp_] := ReplaceAll[exp, {
  FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum + b__FeynCalc`Momentum, Global`mv]] :> 1 / (FeynCalc`SP[a + b, a + b] - Global`mv^2 + I Global`widthv Global`mv),
  FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a_FeynCalc`Momentum - b_FeynCalc`Momentum, Global`mv]] :> 1 / (FeynCalc`SP[a - b, a - b] - Global`mv^2 + I Global`widthv Global`mv),
  FeynCalc`FeynAmpDenominator[FeynCalc`PropagatorDenominator[a__FeynCalc`Momentum, Global`mv]] :> 1 / (FeynCalc`SP[a, a] - Global`mv^2 + I Global`widthv Global`mv)}
];

HazmaSubstituteUnstableMediatorPropagator[exp_] := Module[{},
  If[$HazmaModel === "HazmaScalar", Return[SubstituteUnstableScalarPropagator[exp]]];
  If[$HazmaModel === "HazmaVector", Return[SubstituteUnstableVectorPropagator[exp]]];
  Return[exp]
];

ComputeMandelstamTBound[M_, m1_, m2_, m3_, s_, t_] := Module[{E1, E2},
  E1 = (M^2 + m1^2 - s) / (2 * M);
  E2 = (M^2 + m2^2 - t) / (2 * M);
  Solve[M^2 + 2 * E1 * E2 + m1^2 + m2^2 - m3^2 - 2 * M * (E1 + E2) == 2Sqrt[(E1^2 - m1^2) * (E2^2 - m2^2)], t]
];

KallenLambda[a_, b_, c_] := a^2 + b^2 + c^2 - 2 * (a * b + a * c + b * c);


End[];

EndPackage[]