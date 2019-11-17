(* ::Package:: *)

(* Mathematica Init File *)

Global`$LoadFeynArts = True;
Needs["FeynCalc`"];
$FAVerbose = False;

BeginPackage["HazmaTools`"];

Print["MeV-DM - Scalar Mediator: version 0.2.1"];
Print["by Adam Coogan and Logan A. Morrison"];

Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Utilities.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "States.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Diagrams.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Amplitudes.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "CrossSections.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Widths.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Spectra.m"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName, 2], "Rambo.m"}]];

EndPackage[];
