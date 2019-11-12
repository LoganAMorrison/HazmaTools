(* Mathematica Test File *)

BeginTestSection["TestCrossSections.mt"]

VerificationTest[
  $HazmaModel = "SM";

  SigEEtoGammaGamma = HazmaComputeCrossSection22[{Electron, AntiElectron}, {Photon, Photon}, Q,
    FeynArts`ExcludeParticles -> {Higgs, GoldstoneZBoson, ZBoson}];

  SigEEtoGammaGamma = ReplaceAll[SigEEtoGammaGamma, {FeynCalc`FCGV["ME"] :> FeynCalc`SMP["m_e"],
    FeynCalc`SMP["e"] :> Sqrt[4 * Pi * alphaEM]}];

  SigEEtoGammaGamma = Simplify[SigEEtoGammaGamma, Assumptions -> {Q > 2 SMP["m_e"] > 0}];

  (* peskin page 168 Eqn. 5.106 *)
  SigEEtoGammaGammaPeskin = FullSimplify[
    ReplaceAll[(2 * Pi alphaEM^2) /
        s * (EE / p) * ((EE^2 + p^2 * z^2) / (m^2 + p^2 * (1 - z^2)) + (2 * m^2) /
        (m^2 + p^2 (1 - z^2)) - (2 * m^4) / (m^2 + p^2 * (1 - z^2))^2),
      {s -> Q^2, m -> FeynCalc`SMP["m_e"], EE -> Q / 2, p -> Sqrt[Q^2 / 4 - FeynCalc`SMP["m_e"]^2]}],
    {Q > 0, Q > 2 * FeynCalc`SMP["m_e"], -1 < z < 1}];

  SigEEtoGammaGammaPeskin = (1 / 2) * Integrate[SigEEtoGammaGammaPeskin, {z, -1, 1}, GenerateConditions -> False,
    Assumptions -> {Q > 2 * FeynCalc`SMP["m_e"] > 0, -1 < z < 1}];

  FullSimplify[SigEEtoGammaGammaPeskin - SigEEtoGammaGamma, Assumptions -> {Q > 0, Q > 2 * FeynCalc`SMP["m_e"] > 0}]
  ,
  0
]

EndTestSection[]
