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

GetModelName[] := Module[{},
  If[$HazmaModel === "scalar", Return[{"EFT_MeV_DM_scalar" <> "/" <> "EFT_MeV_DM_scalar"}]];
  If[$HazmaModel === "vector", Return[{"EFT_MeV_DM_vector_no_contact" <> "/" <> "EFT_MeV_DM_vector_no_contact"}]];
  If[$HazmaModel === "SM", Return["SM"]]; (*for testing*)
  If[$HazmaModel === "SMQCD", Return["SMQCD"]]; (*for testing*)
  Throw[$Failed, Abort::InvalidHazmaModel]
];


GetGenericModelName[] := Module[{},
  If[$HazmaModel === "scalar", Return[{"Lorentz", "EFT_MeV_DM_scalar" <> "/" <> "EFT_MeV_DM_scalar"}]];
  If[$HazmaModel === "vector", Return[{"Lorentz", "EFT_MeV_DM_vector_no_contact" <> "/" <> "EFT_MeV_DM_vector_no_contact"}]];
  If[$HazmaModel === "SM", Return[OptionValue[FeynArts`InsertFields, FeynArts`GenericModel]]]; (*for testing*)
  If[$HazmaModel === "SMQCD", Return[OptionValue[FeynArts`InsertFields, FeynArts`GenericModel]]]; (*for testing*)
  Throw[$Failed, Abort::InvalidHazmaModel]
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