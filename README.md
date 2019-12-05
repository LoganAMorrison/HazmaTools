# HazmaTools
`HazmaTools` is a Mathematica package for computing 1->2 widths, 2->2 cross sections and gamma-ray spectra from 2->2 processes.

## Installation
`HazmaTools` requires `FeynCalc` and `FeynArts` patched for use in `FeynCalc` (see https://feyncalc.github.io for details on `FeynCalc`.) To install `HazmaTools`, clone this repo into your Mathematica `Applications` directory. To use the models used in the `Hazma` framework (i.e. the MeV simplified DM models,) copy the directories `HazmaScalar` and `HazmaVector` located in `HazmaTools/Models` into `FeynCalc/FeynArts/Models/`.

To load `HazmaTools` into a Mathematica session, run:

```Mathematica
<<HazmaTools`
```

## Usage

### Setting the model

`HazmaTools` will work with **any** model located inside of your `FeynCalc/FeynArts/Models/` directory.
To set the model, set the global variable `$HazmaModel`. For example, to use models from the `Hazma` framework
(assuming you have copied `HazmaScalar` and `HazmaVector` from `HazmaTools/Models` into `FeynCalc/FeynArts/Models/`)
use:

```mathematica
(* for the simplified MeV DM model with a scalar mediator portal *)
$HazmaTools = "HazmaScalar"
(* for the simplified MeV DM model with a vector mediator portal *)
$HazmaTools = "HazmaVector"
```

As mentioned above, `HazmaTools` will work with any model located in the `FeynCalc/FeynArts/Models/` directory.
I.e., you can use the SM or the MSSM via:

```mathematica
(* use the Standard model: *)
$HazmaTools = "SM"
(* use the Minimal Super Symmetric Model: *)
$HazmaTools = "MSSM"
(* use SM with QCD: *)
$HazmaTools = "SMQCD"
```

Note that switching the model in the middle of a session should work. You may have to run the function twice to get things to work properly. Alternatively, restart the session using `Quit[]` and reload `HazmaTools`.

### Particles
The availible particles for the scalar and vector models are:

```mathematica
DarkMatter/AntiDarkMatter
NeutralPion
ChargedPionM/ChargedPionP
Lepton/AntiLepton (* stands for muon or electron *)
Photon
ScalarMediator (* scalar model only *)
VectorMediator (* vector model only *)
```

These are just aliases for things like `F[1]` for `Neutrino` and so fourth.

The availible particles for `SM` or `SMQCD` are:

```mathematica
(* down-type leptons *)
Electron/AntiElectron
Muon/AntiMuon
Tau/AntiTau

(* neutinos *)
ElectronNeutrino/AntiElectronNeutrino
MuonNeutrino/AntiMuonNeutrino
TauNeutrino/AntiTauNeutrino

(* up-type quarks *)
UpTypeQuark/AntiUpTypeQuark
UpQuark/AntiUpQuark
CharmQuark/AntiCharmQuark
TopQuark/AntiTopQuark

(* down-type quarks *)
DownTypeQuark/AntiDownTypeQuark
DownQuark/AntiDownQuark
StrangeQuark/AntiStrangeQuark
BottomQuark/AntiBottomQuark

(* vector bosons *)
ZBoson
WBosonM/WBosonP
Gluon

(* Higgs and Goldstones *)
Higgs
GoldstoneZBoson
GoldstoneWBosonM/GoldstoneWBosonP

(* ghosts *)
GhostPhoton/AntiGhostPhoton
GhostZBoson/AntiGhostZBoson
GhostWBosonM/AntiGhostWBosonM
GhostWBosonP/AntiGhostWBosonP
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

Alternatively, one can use:

```mathematica
HazmaComputeDNDE[{Lepton, AntiLepton}]
```

which will compute `ScalarMediator -> Lepton + AntiLepton`. Due to factorization of the amplitudes, this will produce the same result as `HazmaComputeDNDE[{DarkMatter, AntiDarkMatter}, {Lepton, AntiLepton},Q]` after the appropriate replacements of the mediator mass with `Q` have been made.
