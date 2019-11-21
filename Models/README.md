# Models
This directory contains the `FeynRules` models used to generate the `FeynArts` model files needed to use HazmaTools.
The folders `ScalarMediator`, `VectorMediator` and so forth contain pre-made model files.

## Building Model Files
To build the model files, open a fresh notebook and run:
```Mathematica
(* Set path to FeynRules. This assumes you have FeynRules installed in standard location, *)
$FeynRulesPath = FileNameJoin[{$UserBaseDirectory, "Applications", "FeynRules"}];

(* To turn off parallelization, use: FR$Parallel = False;*)
<< FeynRules`
```  
Next, load the model files:
```Mathematica
(* This assumes you have HazmaTools installed in standard location. *)
BaseDir = FileNameJoin[{$UserBaseDirectory, "Applications", "HazmaTools", "Models"}];

LoadModel[
 FileNameJoin[{BaseDir, "Info.fr"}],
 FileNameJoin[{BaseDir, "Parameters.fr"}],
 FileNameJoin[{BaseDir, "Particles.fr"}],
 FileNameJoin[{BaseDir, "Tools.fr"}],
 FileNameJoin[{BaseDir, "Lagrangians.fr"}]
]
```
To build the scalar mediator model, use:
```Mathematica
(* Set global flag to turn on scalar mediator model. *)
$IncludeScalarMediator = True;
 
(* Compute the Lagrangian from "Lagrangians.fr" *)
ScalarMediatorLagrangian = Lagrangian;

(* Build FeynArts model files *)
WriteFeynArtsOutput[
  ScalarMediatorLagrangian, 
  MaxParticles -> 5, 
  Output -> NotebookDirectory[] <> "EFT_MeV_DM_scalar"
]
```
One can then build the vector mediator model using:
```Mathematica
(* Turn off scalar mediator model *)
$IncludeScalarMediator = False;

(* Turn on vector mediator model *)
$IncludeVectorMediator = True;

 
(* Compute the Lagrangian from "Lagrangians.fr" *)
VectorMediatorLagrangian = Lagrangian;

(* Build FeynArts model files *)
WriteFeynArtsOutput[
  VectorMediatorLagrangian, 
  MaxParticles -> 5, 
  Output -> NotebookDirectory[] <> "EFT_MeV_DM_vector"
]
```