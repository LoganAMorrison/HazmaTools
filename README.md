# HazmaTools
Mathematica tools for Hazma

To install, clone this repo and navigate to the directory. Then, create a symlink/alias for `HazmaTools.m` to the Mathematica applications directory using:

```code
ln  path/to/HazmaTools.m path/to/Mathematica/Applications
```

Then, to load `HazmaTools` into a Mathematica session, run:

```Mathematica
<<HazmaTools`
```

To set the model, run:

```mathematica
$HazmaTools = "scalar"
```

This will set the model to the scalar-mediator model. Use `"vector"` for the vector mediator model. To squared amplitudes, run:

```mathematica
HazmaCreateAmplitudeSquared[{statex, statexbar}, {statel, statelbar}, IncomingMomenta -> {p1, p2}, OutgoingMomenta -> {p3, p4}]
```

To compute a 2->2 cross section, use:

```mathematica
HazmaComputeCrossSection22[{statex, statexbar}, {statel, statelbar},Q]
```

To compute dNdE, use:

```mathematica
HazmaComputeDNDE[{statex, statexbar}, {statel, statelbar}, Q]
```
