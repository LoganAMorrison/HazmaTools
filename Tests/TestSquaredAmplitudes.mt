(* Mathematica Test File *)

BeginTestSection["TestSquaredAmplitudes.mt"]

VerificationTest[(*1*)
  $HazmaModel = "SM";

  MsqrdBhabha = HazmaComputeAmplitudeSquared[{Electron, AntiElectron}, {Electron, AntiElectron},
    FeynArts`ExcludeParticles -> {Higgs, GoldstoneZBoson, ZBoson},
    FeynCalc`IncomingMomenta -> {p1, p2}, FeynCalc`OutgoingMomenta -> {p3, p4}
  ];
  MsqrdBhabha = ReplaceAll[MsqrdBhabha, {FeynCalc`FCGV["ME"] :> FeynCalc`SMP["m_e"]}];
  FeynCalc`SetMandelstam[s, t, u, p1, p2, -p3, -p4, FeynCalc`SMP["m_e"], FeynCalc`SMP["m_e"],
    FeynCalc`SMP["m_e"], FeynCalc`SMP["m_e"]];

  dSigdzBhabha = (2 * Pi) / (64 * Pi^2 * s) * Simplify[MsqrdBhabha];
  dSigdzBhabha = Normal[Series[dSigdzBhabha, {FeynCalc`SMP["m_e"], 0, 0}]];
  dSigdzBhabha = ReplaceAll[dSigdzBhabha, FeynCalc`SMP["e"] -> Sqrt[4 * Pi * alphaEM]];

  (* peskin pagee 170 problem 5.2 *)
  dSigdzBhabhaPeskin = (Pi * alphaEM^2) / s * (u^2 * (1 / s + 1 / t)^2 + (t / s)^2 + (s / t)^2);
  Simplify[dSigdzBhabha - dSigdzBhabhaPeskin]
  ,
  0
]

VerificationTest[(*2*)
  $HazmaModel = "SM";

  msqrdCompton = HazmaComputeAmplitudeSquared[{Electron, Photon}, {Electron, Photon},
    FeynArts`ExcludeParticles -> {Higgs, GoldstoneZBoson, ZBoson},
    FeynCalc`IncomingMomenta -> {p1, k1},
    FeynCalc`OutgoingMomenta -> {p2, k2}
  ];
  msqrdCompton = ReplaceAll[msqrdCompton, {FeynCalc`FCGV["ME"] :> FeynCalc`SMP["m_e"],
    FeynCalc`SMP["e"] -> Sqrt[4 * Pi * alphaEM]}];

  FeynCalc`SetMandelstam[s, t, u, p1, k1, -p2, -k2, FeynCalc`SMP["m_e"], 0, FeynCalc`SMP["m_e"], 0];
  msqrdCompton = Simplify[msqrdCompton];

  (* Peskin pagee 162 Eqn. 5.87 *)
  msqrdComptonPeskin =
      2 (Sqrt[4 * Pi * alphaEM])^4 * (FeynCalc`SP[p1, k2] / FeynCalc`SP[p1, k1] +
          FeynCalc`SP[p1, k1] / FeynCalc`SP[p1, k2] +
          2 * FeynCalc`SMP["m_e"]^2 * (1 / FeynCalc`SP[p1, k1] - 1 / FeynCalc`SP[p1, k2]) +
          FeynCalc`SMP["m_e"]^4 (1 / FeynCalc`SP[p1, k1] - 1 / FeynCalc`SP[p1, k2])^2);
  FullSimplify[ReplaceAll[
  	msqrdComptonPeskin - msqrdCompton, 
  		{u -> 2 * FeynCalc`SMP["m_e"]^2 - s - t}]]
  ,
  0
]

VerificationTest[(*3*)
  $HazmaModel = "SM";

  msqrdEEtoMuMu = HazmaComputeAmplitudeSquared[{Electron, AntiElectron}, {Muon, AntiMuon},
    FeynArts`ExcludeParticles -> {Higgs, GoldstoneZBoson, ZBoson},
    FeynCalc`IncomingMomenta -> {p1, p2},
    FeynCalc`OutgoingMomenta -> {p3, p4}];

  msqrdEEtoMuMu = ReplaceAll[msqrdEEtoMuMu, {FeynCalc`FCGV["ME"] :> FeynCalc`SMP["m_e"],
    FeynCalc`FCGV["MM"] :> FeynCalc`SMP["m_mu"], FeynCalc`SMP["e"] -> Sqrt[4 * Pi * alphaEM]}];

  SetMandelstam[s, t, u, p1, p2, -p3, -p4, FeynCalc`SMP["m_e"], FeynCalc`SMP["m_e"],
    FeynCalc`SMP["m_mu"], FeynCalc`SMP["m_mu"]];
  msqrdEEtoMuMu = FullSimplify[ReplaceAll[msqrdEEtoMuMu, u -> 2 * FeynCalc`SMP["m_e"]^2 +
      2 * FeynCalc`SMP["m_mu"]^2 - s - t]];

  msqrdEEtoMuMuSchwartz =
      ReplaceAll[(2 * (4 * Pi * alphaEM)^2) / s^2 * (t^2 + u^2 + 4 * s *
          (FeynCalc`SMP["m_e"]^2 + FeynCalc`SMP["m_mu"]^2) -
          2 * (FeynCalc`SMP["m_e"]^2 + FeynCalc`SMP["m_mu"]^2)^2),
        u -> 2 * FeynCalc`SMP["m_e"]^2 + 2 * FeynCalc`SMP["m_mu"]^2 - s - t];

  FullSimplify[msqrdEEtoMuMu - msqrdEEtoMuMuSchwartz]
  ,
  0
]

EndTestSection[]
