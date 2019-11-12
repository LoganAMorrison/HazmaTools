# HazmaTools
Mathematica tools for Hazma

To install, clone this repo into your Mathematica `Applications` directory. Then, to load `HazmaTools` into a Mathematica session, run:

```Mathematica
<<HazmaTools`
```

To set the model, run:

```mathematica
(* set model to scalar mediator model: *)
$HazmaTools = "scalar"
(* set model to vector mediator model: *)
$HazmaTools = "vector"
(* set model to standard model (for testing): *)
$HazmaTools = "SM"
(* set model to standard model QCD (for testing): *)
$HazmaTools = "SMQCD"
```

Note that switching the model in the middle of a session should work. You may have to run the function twice to get things to work properly. Alternatively, restart the session using `Quit[]` and reload `HazmaTools`. The availible particles for the scalar and vector models are: 

```mathematica
DarkMatter/AntiDarkMatter
NeutralPion
ChargedPionM/ChargedPionP
Lepton/AntiLepton (* stands for muon or electron *)
Photon
ScalarMediator (* scalar model only *)
VectorMediator (* vector model only *)
```

To compute an amplitude (or squared amplitude), run:

```mathematica
HazmaComputeAmplitude[{DarkMatter, AntiDarkMatter}, {Lepton, AntiLepton}, IncomingMomenta -> {px, pxbar}, OutgoingMomenta -> {pl, plbar}]
HazmaComputeAmplitudeSquared[{DarkMatter, AntiDarkMatter}, {Lepton, AntiLepton}, IncomingMomenta -> {px, pxbar}, OutgoingMomenta -> {pl, plbar}]
```

To compute a 2->2 cross section, use:

```mathematica
HazmaComputeCrossSection22[{DarkMatter, AntiDarkMatter}, {Lepton, AntiLepton}, Q]
```

To compute a 1->2 width, use:

```mathematica
HazmaComputeWidth[ScalarMediator, {Lepton, AntiLepton}]
```

To compute dNdE, use:

```mathematica
HazmaComputeDNDE[{DarkMatter, AntiDarkMatter}, {Lepton, AntiLepton},Q]
```
