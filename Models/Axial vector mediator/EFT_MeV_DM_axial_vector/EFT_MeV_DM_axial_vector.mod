(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"EFT_MeV_DM_axial_vector",
  Authors -> {"Adam Coogan", "Logan Morrison"},
  Date -> "4/27/2018",
  Institutions -> {"UC Santa Cruz"},
  Emails -> {"acoogan@ucsc.edu", "loanmorr@ucsc.edu"},
  Version -> 1};

FR$ClassesTranslation={};

FR$InteractionOrderPerturbativeExpansion={};

FR$GoldstoneList={};

(*     Declared indices    *)

(*     Declared particles    *)

M$ClassesDescription = {
S[1] == {
    SelfConjugate -> True,
    PropagatorLabel -> ComposedChar["\\pi", Null, "0"],
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> mpi0,
    Indices -> {} },

S[2] == {
    SelfConjugate -> False,
    QuantumNumbers -> {-Q},
    PropagatorLabel -> ComposedChar["\\pi", Null, "-"],
    PropagatorArrow -> Forward,
    PropagatorType -> ScalarDash,
    Mass -> mpi,
    Indices -> {} },

S[3] == {
    SelfConjugate -> False,
    PropagatorLabel -> ComposedChar["K", Null, "0"],
    PropagatorArrow -> Forward,
    PropagatorType -> ScalarDash,
    Mass -> mk0,
    Indices -> {} },

S[4] == {
    SelfConjugate -> False,
    QuantumNumbers -> {-Q},
    PropagatorLabel -> ComposedChar["K", Null, "-"],
    PropagatorArrow -> Forward,
    PropagatorType -> ScalarDash,
    Mass -> mk,
    Indices -> {} },

S[5] == {
    SelfConjugate -> True,
    PropagatorLabel -> "\\eta",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> meta,
    Indices -> {} },

V[1] == {
    SelfConjugate -> True,
    PropagatorLabel -> "\\gamma",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> 0,
    Indices -> {} },

V[2] == {
    SelfConjugate -> True,
    PropagatorLabel -> "Ax",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> mv,
    Indices -> {} },

F[1] == {
    SelfConjugate -> False,
    QuantumNumbers -> {LeptonNumber},
    PropagatorLabel -> "v",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> 0,
    Indices -> {} },

F[2] == {
    SelfConjugate -> False,
    QuantumNumbers -> {-Q, LeptonNumber},
    PropagatorLabel -> "l",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> ml,
    Indices -> {} },

F[3] == {
    SelfConjugate -> False,
    PropagatorLabel -> "\\chi",
    PropagatorArrow -> Forward,
    PropagatorType -> Straight,
    Mass -> mx,
    Indices -> {} }
}


(*        Definitions       *)

GaugeXi[ S[1] ] = 1;
GaugeXi[ S[2] ] = 1;
GaugeXi[ S[3] ] = 1;
GaugeXi[ S[4] ] = 1;
GaugeXi[ S[5] ] = 1;
GaugeXi[ V[1] ] = GaugeXi[A];
GaugeXi[ V[2] ] = GaugeXi[Ax];

mpi0[ ___ ] := mpi0;
mpi[ ___ ] := mpi;
mk0[ ___ ] := mk0;
mk[ ___ ] := mk;
meta[ ___ ] := meta;
mv[ ___ ] := mv;
ml[ ___ ] := ml;
mx[ ___ ] := mx;




(*      Couplings (calculated by FeynRules)      *)

M$CouplingMatrices = {

C[ S[5] , S[5] , S[5] , S[5] ] == {{I*gc1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[5] , S[3] , -S[3] ] == {{(2*I)*b0*gc2*(mdq + 3*msq), 0}, {(6*I)*gc2, 0}, {(-3*I)*gc2, 0}, {(-3*I)*gc2, 0}, {(-3*I)*gc2, 0}, {(-3*I)*gc2, 0}, {(6*I)*gc2, 0}},

C[ S[3] , S[3] , -S[3] , -S[3] ] == {{(-2*I)*b0*gc3*(mdq + msq), 0}, {(2*I)*gc3, 0}, {(-I)*gc3, 0}, {(-I)*gc3, 0}, {(-I)*gc3, 0}, {(-I)*gc3, 0}, {(2*I)*gc3, 0}},

C[ S[3] , -S[3] , S[4] , -S[4] ] == {{I*b0*gc4*(mdq + 2*msq + muq), 0}, {I*gc4, 0}, {I*gc4, 0}, {(-2*I)*gc4, 0}, {(-2*I)*gc4, 0}, {I*gc4, 0}, {I*gc4, 0}},

C[ S[5] , S[5] , S[4] , -S[4] ] == {{(2*I)*b0*gc5*(3*msq + muq), 0}, {(6*I)*gc5, 0}, {(-3*I)*gc5, 0}, {(-3*I)*gc5, 0}, {(-3*I)*gc5, 0}, {(-3*I)*gc5, 0}, {(6*I)*gc5, 0}},

C[ S[4] , S[4] , -S[4] , -S[4] ] == {{(-2*I)*b0*gc6*(msq + muq), 0}, {(2*I)*gc6, 0}, {(-I)*gc6, 0}, {(-I)*gc6, 0}, {(-I)*gc6, 0}, {(-I)*gc6, 0}, {(2*I)*gc6, 0}},

C[ S[5] , S[5] , S[5] , S[1] ] == {{I*gc7, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[3] , -S[3] , S[1] ] == {{(-2*I)*b0*gc8*(mdq - msq), 0}, {(3*I)*gc8, 0}, {(3*I)*gc8, 0}, {(-6*I)*gc8, 0}, {(-6*I)*gc8, 0}, {(3*I)*gc8, 0}, {(3*I)*gc8, 0}},

C[ S[5] , S[4] , -S[4] , S[1] ] == {{(2*I)*b0*gc9*(msq - muq), 0}, {(3*I)*gc9, 0}, {(3*I)*gc9, 0}, {(-6*I)*gc9, 0}, {(-6*I)*gc9, 0}, {(3*I)*gc9, 0}, {(3*I)*gc9, 0}},

C[ S[5] , S[5] , S[1] , S[1] ] == {{I*gc10, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[1] , S[1] ] == {{(2*I)*b0*gc11*(3*mdq + msq), 0}, {(2*I)*gc11, 0}, {(-I)*gc11, 0}, {(-I)*gc11, 0}, {(-I)*gc11, 0}, {(-I)*gc11, 0}, {(2*I)*gc11, 0}},

C[ S[4] , -S[4] , S[1] , S[1] ] == {{(2*I)*b0*gc12*(msq + 3*muq), 0}, {(2*I)*gc12, 0}, {(-I)*gc12, 0}, {(-I)*gc12, 0}, {(-I)*gc12, 0}, {(-I)*gc12, 0}, {(2*I)*gc12, 0}},

C[ S[5] , S[1] , S[1] , S[1] ] == {{I*gc13, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[1] , S[1] , S[1] , S[1] ] == {{I*gc14, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , -S[3] , -S[4] , S[2] ] == {{(-I)*b0*gc15*(mdq - 2*msq + muq), 0}, {(3*I)*gc15, 0}, {(3*I)*gc15, 0}, {(-6*I)*gc15, 0}, {(-6*I)*gc15, 0}, {(3*I)*gc15, 0}, {(3*I)*gc15, 0}},

C[ -S[3] , -S[4] , S[1] , S[2] ] == {{I*b0*gc16*(mdq - muq), 0}, {0, 0}, {(3*I)*gc16, 0}, {(-3*I)*gc16, 0}, {(-3*I)*gc16, 0}, {(3*I)*gc16, 0}, {0, 0}},

C[ S[5] , S[3] , S[4] , -S[2] ] == {{(-I)*b0*gc17*(mdq - 2*msq + muq), 0}, {(3*I)*gc17, 0}, {(3*I)*gc17, 0}, {(-6*I)*gc17, 0}, {(-6*I)*gc17, 0}, {(3*I)*gc17, 0}, {(3*I)*gc17, 0}},

C[ S[3] , S[4] , S[1] , -S[2] ] == {{I*b0*gc18*(mdq - muq), 0}, {0, 0}, {(3*I)*gc18, 0}, {(-3*I)*gc18, 0}, {(-3*I)*gc18, 0}, {(3*I)*gc18, 0}, {0, 0}},

C[ S[5] , S[5] , S[2] , -S[2] ] == {{I*gc19, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[2] , -S[2] ] == {{I*b0*gc20*(2*mdq + msq + muq), 0}, {I*gc20, 0}, {(-2*I)*gc20, 0}, {I*gc20, 0}, {I*gc20, 0}, {(-2*I)*gc20, 0}, {I*gc20, 0}},

C[ S[4] , -S[4] , S[2] , -S[2] ] == {{I*b0*gc21*(mdq + msq + 2*muq), 0}, {I*gc21, 0}, {(-2*I)*gc21, 0}, {I*gc21, 0}, {I*gc21, 0}, {(-2*I)*gc21, 0}, {I*gc21, 0}},

C[ S[5] , S[1] , S[2] , -S[2] ] == {{I*gc22, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[1] , S[1] , S[2] , -S[2] ] == {{I*b0*gc23*(mdq + muq), 0}, {(2*I)*gc23, 0}, {(-I)*gc23, 0}, {(-I)*gc23, 0}, {(-I)*gc23, 0}, {(-I)*gc23, 0}, {(2*I)*gc23, 0}},

C[ S[2] , S[2] , -S[2] , -S[2] ] == {{(-2*I)*b0*gc24*(mdq + muq), 0}, {(2*I)*gc24, 0}, {(-I)*gc24, 0}, {(-I)*gc24, 0}, {(-I)*gc24, 0}, {(-I)*gc24, 0}, {(2*I)*gc24, 0}},

C[ S[4] , -S[4] , V[1] , V[1] ] == {{(2*I)*qe^2, 0}},

C[ S[2] , -S[2] , V[1] , V[1] ] == {{(2*I)*qe^2, 0}},

C[ S[3] , -S[3] , V[2] , V[2] ] == {{(-2*I)*(gadd - gass)^2, 0}},

C[ S[4] , -S[4] , V[2] , V[2] ] == {{(-2*I)*(gass - gauu)^2, 0}},

C[ S[2] , -S[2] , V[2] , V[2] ] == {{(-2*I)*(gadd - gauu)^2, 0}},

C[ S[5] , S[3] , -S[3] , V[2] ] == {{-2*gc30, 0}, {gc30, 0}, {gc30, 0}},

C[ S[5] , S[4] , -S[4] , V[2] ] == {{-2*gc31, 0}, {gc31, 0}, {gc31, 0}},

C[ S[3] , -S[3] , S[1] , V[2] ] == {{-gc32, 0}, {-gc32, 0}, {2*gc32, 0}},

C[ S[3] , S[4] , -S[2] , V[2] ] == {{(gadd + gass - 2*gauu)*gc33, 0}, {(-2*gadd + gass + gauu)*gc33, 0}, {(gadd - 2*gass + gauu)*gc33, 0}},

C[ -S[3] , -S[4] , S[2] , V[2] ] == {{(gadd + gass - 2*gauu)*gc34, 0}, {(-2*gadd + gass + gauu)*gc34, 0}, {(gadd - 2*gass + gauu)*gc34, 0}},

C[ S[4] , -S[4] , V[1] ] == {{(-I)*gc35, 0}, {I*gc35, 0}},

C[ S[4] , -S[4] , S[1] , V[2] ] == {{-gc36, 0}, {-gc36, 0}, {2*gc36, 0}},

C[ S[1] , S[2] , -S[2] , V[2] ] == {{-2*gc37, 0}, {gc37, 0}, {gc37, 0}},

C[ S[2] , -S[2] , V[1] ] == {{(-I)*gc38, 0}, {I*gc38, 0}},

C[ -F[2] , F[2] , V[1] ] == {{I*gc39, 0}, {I*gc39, 0}},

C[ -F[2] , F[2] , V[2] ] == {{I*gc40L, 0}, {I*gc40R, 0}},

C[ -F[3] , F[3] , V[2] ] == {{I*gc41L, 0}, {I*gc41R, 0}}

}

(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)

(* Parameter replacement lists (These lists were created by FeynRules) *)

(* FA Couplings *)

M$FACouplings = {
     gc1 -> (b0*(mdq + 16*msq + muq))/(9*fpi^2),
     gc2 -> 1/(12*fpi^2),
     gc3 -> -1/(3*fpi^2),
     gc4 -> 1/(6*fpi^2),
     gc5 -> 1/(12*fpi^2),
     gc6 -> -1/(3*fpi^2),
     gc7 -> -(b0*(mdq - muq))/(3*Sqrt[3]*fpi^2),
     gc8 -> 1/(12*Sqrt[3]*fpi^2),
     gc9 -> -1/(12*Sqrt[3]*fpi^2),
     gc10 -> (b0*(mdq + muq))/(3*fpi^2),
     gc11 -> 1/(12*fpi^2),
     gc12 -> 1/(12*fpi^2),
     gc13 -> -((b0*(mdq - muq))/(Sqrt[3]*fpi^2)),
     gc14 -> (b0*(mdq + muq))/fpi^2,
     gc15 -> -1/(6*Sqrt[6]*fpi^2),
     gc16 -> -1/(6*Sqrt[2]*fpi^2),
     gc17 -> -1/(6*Sqrt[6]*fpi^2),
     gc18 -> -1/(6*Sqrt[2]*fpi^2),
     gc19 -> (b0*(mdq + muq))/(3*fpi^2),
     gc20 -> 1/(6*fpi^2),
     gc21 -> 1/(6*fpi^2),
     gc22 -> -(b0*(mdq - muq))/(3*Sqrt[3]*fpi^2),
     gc23 -> 1/(3*fpi^2),
     gc24 -> -1/(3*fpi^2),
     gc30 -> (-gadd + gass)/(Sqrt[3]*fpi),
     gc31 -> (gass - gauu)/(Sqrt[3]*fpi),
     gc32 -> (-gadd + gass)/(3*fpi),
     gc33 -> Sqrt[2]/(3*fpi),
     gc34 -> Sqrt[2]/(3*fpi),
     gc35 -> qe,
     gc36 -> (-gass + gauu)/(3*fpi),
     gc37 -> (2*(gadd - gauu))/(3*fpi),
     gc38 -> qe,
     gc39 -> -qe,
     gc40L -> -gall,
     gc40R -> gall,
     gc41L -> -gaxx,
     gc41R -> gaxx};
