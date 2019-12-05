#!/usr/bin/env wolframscript

(* Wolfram Script *)
(* :Author: loganmorrison *)
(* :Date: 2019-12-03 *)

$FeynRulesPath = SetDirectory[FileNameJoin[{$UserBaseDirectory, "Applications", "FeynRules"}]];
<< FeynRules`;

ModelDirectory = FileNameJoin[{$UserBaseDirectory, "Applications", "HazmaTools", "Models"}];
LoadModel[
  FileNameJoin[{ModelDirectory, "Info.fr"}],
  FileNameJoin[{ModelDirectory, "Tools.fr"}],
  FileNameJoin[{ModelDirectory, "Parameters.fr"}],
  FileNameJoin[{ModelDirectory, "Particles.fr"}],
  FileNameJoin[{ModelDirectory, "Lagrangians.fr"}]
];

$IncludeScalarMediator = True;
SetDirectory[FileNameJoin[{$UserBaseDirectory, "Applications", "FeynCalc", "FeynArts", "Models"}]];
WriteFeynArtsOutput[Lagrangian, MaxParticles -> 5, Output -> "HazmaScalar", DiracIndices -> False, CouplingRename -> False];
$IncludeScalarMediator = False;

$IncludeVectorMediator = True;
SetDirectory[FileNameJoin[{$UserBaseDirectory, "Applications", "FeynCalc", "FeynArts", "Models"}]];
WriteFeynArtsOutput[Lagrangian, MaxParticles -> 5, Output -> "HazmaVector", DiracIndices -> False, CouplingRename -> False];
$IncludeVectorMediator = False;