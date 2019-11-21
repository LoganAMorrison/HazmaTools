(* Mathematica Package *)

(* :Title: Rambo *)
(* :Context: HazmaTools` *)
(* :Author: Adam Coogan and Logan Morrison *)
(* :Date: 2019-11-16 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 Adam Coogan and Logan Morrison *)
(* :Keywords:  Quantum Field Theory *)
(* :Discussion: File for integrating phase space *)

BeginPackage["HazmaTools`"];
(* Exported symbols added here with SymbolName::usage *)

RamboGeneratePhaseSpace::usage = "Generate phase space points given center-of-mass energy `cme` and final state \
particle masses `masses`.";

RamboIntegratePhaseSpace::usage = "Integrate over phase space given center-of-mass energy `cme` and final state \
particle masses `masses`.";

RamboCrossSection::usage = "Compute cross-section given center-of-mass energy `cme`, initial-state-particle masses \
`ispMasses`, final-state-particle masses `fspMasses`.";

RamboDecayWidth::usage = "Compute decay width given parent particle mass `M` and final-state-particle masses
`fspMasses`.";

Parallel::usage = "Option for Rambo functions to use parallelization.";
NumPoints::usage = "Option for Rambo functions to specify number of points to generate.";

Begin["`Private`"];

(*
  Core Functions
  --------------
  -`Metric`: 4D west-coast Minkowski metric
  -`MassSquared`: computes scalar product of four-momenta with itself
  -`GenerateRandomFourMomentum`: generate isotropic massless four-momenta
  -`BoostSingleFourMomentum`: boost a four-vector by a given boost vector
  -`BoostFourMomenta`: boost all four momenta into center of mass frame
  -`ComputeScaleFactor`: compute factor to correct masses of four momenta
  -`CorrectMasses`: correct four-vectors to have the correct mass
  -`GeneratePhaseSpacePoint` generate a set of four-momenta in center of mass frame with correct masses
*)

SetAttributes[MassSquared, NumericFunction];
SetAttributes[GenerateRandomFourMomentum, NumericFunction];
SetAttributes[BoostSingleFourMomentum, NumericFunction];
SetAttributes[BoostFourMomenta, NumericFunction];
SetAttributes[ComputeScaleFactor, NumericFunction];
SetAttributes[CorrectMasses, NumericFunction];
SetAttributes[GeneratePhaseSpacePoint, NumericFunction];

Metric = DiagonalMatrix[{1, -1, -1, -1}];

MassSquared[fm_?VectorQ /; Length[fm] == 4] := Dot[fm, Metric, fm];

GenerateRandomFourMomentum[] := Module[{e, c, phi},
  c = 2.0 * RandomReal[] - 1.0;
  phi = 2 * Pi * RandomReal[];
  e = -Log[RandomReal[] * RandomReal[]];
  {e, e * Sqrt[1 - c^2] * Cos[phi], e * Sqrt[1 - c^2] * Sin[phi], e * c}
];

BoostSingleFourMomentum[fm_?VectorQ, x_?NumericQ, b_?VectorQ, gamma_?NumericQ, a_?NumericQ] := Module[{dotProd},
  dotProd = b[[1]] * fm[[2]] + b[[2]] * fm[[3]] + b[[3]] * fm[[4]];
  {
    x * (gamma * fm[[1]] + dotProd),
    x * (fm[[2]] + b[[1]] * fm[[1]] + a * dotProd * b[[1]]),
    x * (fm[[3]] + b[[2]] * fm[[1]] + a * dotProd * b[[2]]),
    x * (fm[[4]] + b[[3]] * fm[[1]] + a * dotProd * b[[3]])
  }
];

BoostFourMomenta[fms_?MatrixQ, cme_?NumericQ] := Module[{sumQs, massQ, b, x, gamma, a, n},
  sumQs = Total[fms];
  massQ = Re[Sqrt[MassSquared[sumQs]]];

  b = {-sumQs[[2]] / massQ, -sumQs[[3]] / massQ, -sumQs[[4]] / massQ};

  x = cme / massQ;
  gamma = sumQs[[1]] / massQ;
  a = 1 / (1 + gamma);
  n = Length[fms];

  (* return list of boosted four-momenta and weight *)
  {
    Table[BoostSingleFourMomentum[fm, x, b, gamma, a], {fm, fms}],
    (Pi / 2)^(n - 1) * cme^(2 * n - 4) / Factorial[n - 2] / Factorial[n - 1] * (2 * Pi)^(4 - 3 * n)
  }
];

ComputeScaleFactor[fms_?MatrixQ, cme_?NumericQ, masses_?VectorQ] := Module[{xi, dxi, isDone, tol, maxIter, f, df, delf, iter},
  xi = Sqrt[1 - Total[masses]^2 / cme^2];

  tol = 10^-4;
  maxIter = 50;
  isDone = False;
  iter = 0;
  While[isDone == False,
    {f, df} = {-cme, 0.0};
    For[i = 1, i <= Length[masses], i++,
      delf = Sqrt[masses[[i]]^2 + xi^2 * fms[[i, 1]]^2];
      f += delf;
      df += xi * fms[[i, 1]]^2 / delf;
    ];
    dxi = -f / df;
    xi += dxi;

    iter += 1;
    If[Or[Abs[dxi] < tol, iter > maxIter], isDone = True]
  ];
  xi
];

CorrectMasses[fms_?MatrixQ, weight_?NumericQ, cme_?NumericQ, masses_?VectorQ] := Module[{t1, t2, t3, xi, newFM, n},

  xi = ComputeScaleFactor[fms, cme, masses];
  newFM = Table[
    {
      Sqrt[masses[[i]]^2 + xi^2 * fms[[i, 1]]^2],
      xi * fms[[i, 2]],
      xi * fms[[i, 3]],
      xi * fms[[i, 4]]
    },
    {i, Length[masses]}
  ];

  t1 = Sum[Sqrt[fm[[2]]^2 + fm[[3]]^2 + fm[[4]]^2] / cme, {fm, newFM}];
  t2 = Sum[(fm[[2]]^2 + fm[[3]]^2 + fm[[4]]^2) / fm[[1]], {fm, newFM}];
  t3 = Product[Sqrt[fm[[2]]^2 + fm[[3]]^2 + fm[[4]]^2] / fm[[1]], {fm, newFM}];

  n = Length[masses];
  t1 = t1^(2 * n - 3);
  t2 = 1 / t2;

  {
    newFM,
    weight * t1 * t2 * t3 * cme
  }

];

GeneratePhaseSpacePoint[cme_?NumericQ, masses_?VectorQ] := Module[{fms, weight},
  fms = Table[GenerateRandomFourMomentum[], {Length[masses]}];
  {fms, weight} = BoostFourMomenta[fms, cme];
  {fms, weight} = CorrectMasses[fms, weight, cme, masses];
  {fms, weight}
];

GeneratePhaseSpacePoint[cme_?NumericQ, masses_?VectorQ, msqrd_] := Module[{fms, weight},
  fms = Table[GenerateRandomFourMomentum[], {Length[masses]}];
  {fms, weight} = BoostFourMomenta[fms, cme];
  {fms, weight} = CorrectMasses[fms, weight, cme, masses];
  {fms, weight * msqrd[fms]}
];

(*
  API Functions
  -------------
  -`RamboGeneratePhaseSpace`: Generate many phase space points
  -`RamboIntegratePhaseSpace`: Integrate over a phase space
  -`RamboCrossSection`: Compute a cross section for 2 -> N
  -`RamboDecayWidth`: Compute a width for 1 -> N
*)

SetAttributes[RamboGeneratePhaseSpace, NumericFunction];
SetAttributes[RamboIntegratePhaseSpace, NumericFunction];
SetAttributes[RamboCrossSection, NumericFunction];
SetAttributes[RamboDecayWidth, NumericFunction];

Options[RamboGeneratePhaseSpace] := {
  NumPoints -> 10000,
  Parallel -> True
};

RamboGeneratePhaseSpace[cme_?NumericQ, masses_?VectorQ, OptionsPattern[]] := Module[{},
  If[OptionValue[Parallel],
    ParallelTable[GeneratePhaseSpacePoint[cme, masses], OptionValue[NumPoints]],
    Table[GeneratePhaseSpacePoint[cme, masses], OptionValue[NumPoints]]
  ]
];

RamboGeneratePhaseSpace[cme_?NumericQ, masses_?VectorQ, msqrd_, OptionsPattern[]] := Module[{},
  If[OptionValue[Parallel],
    ParallelTable[GeneratePhaseSpacePoint[cme, masses, msqrd], OptionValue[NumPoints]],
    Table[GeneratePhaseSpacePoint[cme, masses, msqrd], OptionValue[NumPoints]]
  ]
];

RamboIntegratePhaseSpace[cme_?NumericQ, masses_?VectorQ, OptionsPattern[RamboGeneratePhaseSpace]] := Module[{points, weights},
  points = RamboGeneratePhaseSpace[cme, masses,
    Parallel -> OptionValue[Parallel],
    NumPoints -> OptionValue[NumPoints]];
  weights = If[OptionValue[Parallel],
    ParallelTable[point[[2]], {point, points}],
    Table[point[[2]], {point, points}]
  ];
  {Mean[weights], Sqrt[Variance[weights] / OptionValue[NumPoints]]}
];

RamboIntegratePhaseSpace[cme_?NumericQ, masses_?VectorQ, msqrd_, OptionsPattern[RamboGeneratePhaseSpace]] := Module[{points, weights},
  points = RamboGeneratePhaseSpace[cme, masses, msqrd,
    Parallel -> OptionValue[Parallel],
    NumPoints -> OptionValue[NumPoints]
  ];
  weights = If[OptionValue[Parallel],
    ParallelTable[point[[2]], {point, points}],
    Table[point[[2]], {point, points}]
  ];
  {Mean[weights], Sqrt[Variance[weights] / OptionValue[NumPoints]]}
];

RamboCrossSection[cme_?NumericQ, ispMasses_?VectorQ, fspMasses_?VectorQ, OptionsPattern[RamboGeneratePhaseSpace]] := Module[{integral, err, m1, m2, p, pf},
  {integral, err} = RamboIntegratePhaseSpace[cme, fspMasses,
    Parallel -> OptionValue[Parallel],
    NumPoints -> OptionValue[NumPoints]
  ];
  m1 = ispMasses[[1]];
  m2 = ispMasses[[2]];

  p = Sqrt[(m1 - m2 - cme) * (m1 + m2 - cme) * (m1 - m2 + cme) * (m1 + m2 + cme)] / (2 * cme);

  pf = 1 / (4 * p * cme);

  {integral * pf, err * pf}
];

RamboCrossSection[cme_?NumericQ, ispMasses_?VectorQ, fspMasses_?VectorQ, msqrd_, OptionsPattern[RamboGeneratePhaseSpace]] := Module[{integral, err, m1, m2, p, pf},
  {integral, err} = RamboIntegratePhaseSpace[cme, fspMasses, msqrd,
    Parallel -> OptionValue[Parallel],
    NumPoints -> OptionValue[NumPoints]
  ];
  m1 = ispMasses[[1]];
  m2 = ispMasses[[2]];

  p = Sqrt[(m1 - m2 - cme) * (m1 + m2 - cme) * (m1 - m2 + cme) * (m1 + m2 + cme)] / (2 * cme);

  pf = 1 / (4 * p * cme);

  {integral * pf, err * pf}
];

RamboDecayWidth[M_?NumericQ, fspMasses_?VectorQ, OptionsPattern[RamboGeneratePhaseSpace]] := Module[{integral, err},
  {integral, err} = RamboIntegratePhaseSpace[M, fspMasses,
    Parallel -> OptionValue[Parallel],
    NumPoints -> OptionValue[NumPoints]
  ];
  {integral / (2 * M), err / (2 * M)}
];

RamboDecayWidth[M_?NumericQ, fspMasses_?VectorQ, msqrd_, OptionsPattern[RamboGeneratePhaseSpace]] := Module[{integral, err},
  {integral, err} = RamboIntegratePhaseSpace[M, fspMasses, msqrd,
    Parallel -> OptionValue[Parallel],
    NumPoints -> OptionValue[NumPoints]
  ];
  {integral / (2 * M), err / (2 * M)}
];

End[]; (* `Private` *)

EndPackage[]