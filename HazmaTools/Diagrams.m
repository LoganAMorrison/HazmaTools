(* Mathematica Package *)

(* :Title: Diagrams *)
(* :Context: HazmaTools` *)
(* :Author: Adam Coogan and Logan Morrison *)
(* :Date: 2019-11-11 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 12.0 *)
(* :Copyright: (c) 2019 Adam Coogan and Logan Morrison *)
(* :Keywords: Hazma, Quantum Field Theory, FeynCalc, FeynArts *)
(* :Discussion: File for computing diagrams from FeynArts in the Hazma framework *)

BeginPackage["HazmaTools`"];

HazmaGenerateDiagrams::usage = "HazmaGenerateDiagrams[inStates, outStates] Compute the diagrams\
associated with the process inStates->outStates.";

Begin["`Private`"];

GetModelName::InvalidModel = "Invalid model: `1`.";

GetModelName[] := Module[{dir, type, model},
  model = None;
  dir = Directory[];
  SetDirectory[FileNameJoin[{$UserBaseDirectory, "Applications", "FeynCalc", "FeynArts", "Models"}]];
  (* If model is a directory, need to use $HazmaModel/$HazmaModel *)
  type = FileType[$HazmaModel];
  If[SameQ[type, Directory], model = {FileNameJoin[{$HazmaModel, $HazmaModel}]}];
  (* If there is a bare .mod file, we just use defaults *)
  type = FileType[$HazmaModel <> ".mod"];
  If[SameQ[type, File], model = $HazmaModel];
  (* If gmodel is None, no model was found *)
  If[UnsameQ[model, None],
    SetDirectory[dir]; Return[model],
    Message[GetModelName::InvalidModel, $HazmaModel]; $Failed
  ]
];

GetGenericModelName[] := Module[{dir, type, gmodel},
  gmodel = None;
  dir = Directory[];
  SetDirectory[FileNameJoin[{$UserBaseDirectory, "Applications", "FeynCalc", "FeynArts", "Models"}]];
  (* If model is a directory, need to add Lorentz and $HazmaModel/$HazmaModel *)
  type = FileType[$HazmaModel];
  If[SameQ[type, Directory], gmodel = {"Lorentz", FileNameJoin[{$HazmaModel, $HazmaModel}]}];
  (* If there is a bare .mod file, we just use defaults *)
  type = FileType[$HazmaModel <> ".mod"];
  If[SameQ[type, File], gmodel = Return[OptionValue[FeynArts`InsertFields, FeynArts`GenericModel]]];
  (* If gmodel is None, no model was found *)
  If[UnsameQ[gmodel, None],
    SetDirectory[dir]; Return[gmodel],
    Message[GetModelName::InvalidModel, $HazmaModel]; $Failed
  ]
];

Options[HazmaGenerateDiagrams] = {
  FeynArts`LoopNumber -> 0,
  FeynArts`Adjacencies -> {3, 4, 5},
  FeynArts`ExcludeParticles -> {},
  FeynArts`ExcludeTopologies -> {FeynArts`Tadpoles, FeynArts`SelfEnergies, FeynArts`WFCorrections},
  FeynArts`InsertionLevel -> {FeynArts`Particles},
  FeynArts`Paint -> False,
  FeynArts`ColumnsXRows -> OptionValue[FeynArts`Paint, FeynArts`ColumnsXRows]
};

HazmaGenerateDiagrams[inStates_, outStates_, OptionsPattern[]] := Block[{tops, diags},
  tops = FeynArts`CreateTopologies[
    OptionValue[FeynArts`LoopNumber],
    Length[inStates] -> Length[outStates],
    FeynArts`ExcludeTopologies -> OptionValue[FeynArts`ExcludeTopologies],
    FeynArts`Adjacencies -> OptionValue[FeynArts`Adjacencies]
  ];
  diags = FeynArts`InsertFields[
    tops,
    inStates -> outStates,
    FeynArts`InsertionLevel -> OptionValue[FeynArts`InsertionLevel],
    FeynArts`Model -> GetModelName[],
    FeynArts`GenericModel -> GetGenericModelName[],
    FeynArts`ExcludeParticles -> OptionValue[FeynArts`ExcludeParticles]
  ];
  If[OptionValue[FeynArts`Paint], FeynArts`Paint[diags, FeynArts`ColumnsXRows -> OptionValue[FeynArts`ColumnsXRows]], FeynArts`SheetHeader -> None];
  diags
];

End[];

EndPackage[]