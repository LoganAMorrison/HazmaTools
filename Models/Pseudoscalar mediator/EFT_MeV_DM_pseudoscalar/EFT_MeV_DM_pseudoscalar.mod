(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"EFT_subGeV_DM_P",
  Authors -> {"Adam Coogan", "Logan Morrison"},
  Date -> "7/10/2018",
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

S[7] == {
    SelfConjugate -> True,
    PropagatorLabel -> p,
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> mp,
    Indices -> {} },

V[1] == {
    SelfConjugate -> True,
    PropagatorLabel -> A,
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> 0,
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

GaugeXi[ V[1] ] = GaugeXi[A];

mpi0[ ___ ] := mpi0;
mpi[ ___ ] := mpi;
mk0[ ___ ] := mk0;
mk[ ___ ] := mk;
meta[ ___ ] := meta;
mp[ ___ ] := mp;
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

C[ S[5] , S[5] , S[5] , S[7] ] == {{I*gc7, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[3] , -S[3] , S[7] ] == {{(-2*I)*b0*gc8*((mdq - msq)*sinbeta*vh + 6*cosbeta*fpi*(gpGG*msq + gpss*vh)), 0}, {(3*I)*gc8*sinbeta*vh, 0}, {(3*I)*gc8*sinbeta*vh, 0}, {(-6*I)*gc8*sinbeta*vh, 0}, {(-6*I)*gc8*sinbeta*vh, 0}, {(3*I)*gc8*sinbeta*vh, 0}, {(3*I)*gc8*sinbeta*vh, 0}},

C[ S[5] , S[4] , -S[4] , S[7] ] == {{(2*I)*b0*gc9*((msq - muq)*sinbeta*vh + 6*cosbeta*fpi*(gpGG*msq + gpss*vh)), 0}, {(3*I)*gc9*sinbeta*vh, 0}, {(3*I)*gc9*sinbeta*vh, 0}, {(-6*I)*gc9*sinbeta*vh, 0}, {(-6*I)*gc9*sinbeta*vh, 0}, {(3*I)*gc9*sinbeta*vh, 0}, {(3*I)*gc9*sinbeta*vh, 0}},

C[ S[5] , S[5] , S[5] , S[1] ] == {{I*gc10, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[3] , -S[3] , S[1] ] == {{(2*I)*b0*gc11*(cosbeta*(-mdq + msq)*vh + 6*fpi*sinbeta*(gpGG*msq + gpss*vh)), 0}, {(3*I)*cosbeta*gc11*vh, 0}, {(3*I)*cosbeta*gc11*vh, 0}, {(-6*I)*cosbeta*gc11*vh, 0}, {(-6*I)*cosbeta*gc11*vh, 0}, {(3*I)*cosbeta*gc11*vh, 0}, {(3*I)*cosbeta*gc11*vh, 0}},

C[ S[5] , S[4] , -S[4] , S[1] ] == {{(-2*I)*b0*gc12*(cosbeta*(-msq + muq)*vh + 6*fpi*sinbeta*(gpGG*msq + gpss*vh)), 0}, {(3*I)*cosbeta*gc12*vh, 0}, {(3*I)*cosbeta*gc12*vh, 0}, {(-6*I)*cosbeta*gc12*vh, 0}, {(-6*I)*cosbeta*gc12*vh, 0}, {(3*I)*cosbeta*gc12*vh, 0}, {(3*I)*cosbeta*gc12*vh, 0}},

C[ S[5] , S[5] , S[7] , S[1] ] == {{I*gc13, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[7] , S[1] ] == {{(-2*I)*b0*gc14*(-(cosbeta*(3*mdq + msq)*sinbeta*vh) + 2*cosbeta^2*fpi*(gpGG*(2*mdq + msq) + (2*gpdd + gpss)*vh) - 2*fpi*sinbeta^2*(gpGG*(2*mdq + msq) + (2*gpdd + gpss)*vh)), 0}, {(2*I)*cosbeta*gc14*sinbeta*vh, 0}, {(-I)*cosbeta*gc14*sinbeta*vh, 0}, {(-I)*cosbeta*gc14*sinbeta*vh, 0}, {(-I)*cosbeta*gc14*sinbeta*vh, 0}, {(-I)*cosbeta*gc14*sinbeta*vh, 0}, {(2*I)*cosbeta*gc14*sinbeta*vh, 0}},

C[ S[4] , -S[4] , S[7] , S[1] ] == {{(2*I)*b0*gc15*(cosbeta*(msq + 3*muq)*sinbeta*vh + 2*cosbeta^2*fpi*(gpGG*(msq + 2*muq) + (gpss + 2*gpuu)*vh) - 2*fpi*sinbeta^2*(gpGG*(msq + 2*muq) + (gpss + 2*gpuu)*vh)), 0}, {(2*I)*cosbeta*gc15*sinbeta*vh, 0}, {(-I)*cosbeta*gc15*sinbeta*vh, 0}, {(-I)*cosbeta*gc15*sinbeta*vh, 0}, {(-I)*cosbeta*gc15*sinbeta*vh, 0}, {(-I)*cosbeta*gc15*sinbeta*vh, 0}, {(2*I)*cosbeta*gc15*sinbeta*vh, 0}},

C[ S[5] , S[5] , S[1] , S[1] ] == {{I*gc16, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[1] , S[1] ] == {{(2*I)*b0*gc17*(cosbeta*(3*mdq + msq)*vh + 4*fpi*sinbeta*(gpGG*(2*mdq + msq) + (2*gpdd + gpss)*vh)), 0}, {(2*I)*cosbeta*gc17*vh, 0}, {(-I)*cosbeta*gc17*vh, 0}, {(-I)*cosbeta*gc17*vh, 0}, {(-I)*cosbeta*gc17*vh, 0}, {(-I)*cosbeta*gc17*vh, 0}, {(2*I)*cosbeta*gc17*vh, 0}},

C[ S[4] , -S[4] , S[1] , S[1] ] == {{(2*I)*b0*gc18*(cosbeta*(msq + 3*muq)*vh - 4*fpi*sinbeta*(gpGG*(msq + 2*muq) + (gpss + 2*gpuu)*vh)), 0}, {(2*I)*cosbeta*gc18*vh, 0}, {(-I)*cosbeta*gc18*vh, 0}, {(-I)*cosbeta*gc18*vh, 0}, {(-I)*cosbeta*gc18*vh, 0}, {(-I)*cosbeta*gc18*vh, 0}, {(2*I)*cosbeta*gc18*vh, 0}},

C[ S[5] , S[7] , S[1] , S[1] ] == {{I*gc19, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[1] , S[1] , S[1] ] == {{I*gc20, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[7] , S[1] , S[1] , S[1] ] == {{I*gc21, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[1] , S[1] , S[1] , S[1] ] == {{I*gc22, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , -S[3] , -S[4] , S[2] ] == {{(-I)*b0*gc23*(mdq - 2*msq + muq), 0}, {(3*I)*gc23, 0}, {(3*I)*gc23, 0}, {(-6*I)*gc23, 0}, {(-6*I)*gc23, 0}, {(3*I)*gc23, 0}, {(3*I)*gc23, 0}},

C[ -S[3] , -S[4] , S[7] , S[2] ] == {{(-I)*b0*gc24*((-mdq + muq)*sinbeta*vh + 4*cosbeta*fpi*(gpGG*(mdq + msq + muq) + (gpdd + gpss + gpuu)*vh)), 0}, {0, 0}, {(3*I)*gc24*sinbeta*vh, 0}, {(-3*I)*gc24*sinbeta*vh, 0}, {(-3*I)*gc24*sinbeta*vh, 0}, {(3*I)*gc24*sinbeta*vh, 0}, {0, 0}},

C[ -S[3] , -S[4] , S[1] , S[2] ] == {{I*b0*gc25*(cosbeta*(mdq - muq)*vh + 4*fpi*sinbeta*(gpGG*(mdq + msq + muq) + (gpdd + gpss + gpuu)*vh)), 0}, {0, 0}, {(3*I)*cosbeta*gc25*vh, 0}, {(-3*I)*cosbeta*gc25*vh, 0}, {(-3*I)*cosbeta*gc25*vh, 0}, {(3*I)*cosbeta*gc25*vh, 0}, {0, 0}},

C[ S[5] , S[3] , S[4] , -S[2] ] == {{(-I)*b0*gc26*(mdq - 2*msq + muq), 0}, {(3*I)*gc26, 0}, {(3*I)*gc26, 0}, {(-6*I)*gc26, 0}, {(-6*I)*gc26, 0}, {(3*I)*gc26, 0}, {(3*I)*gc26, 0}},

C[ S[3] , S[4] , S[7] , -S[2] ] == {{(-I)*b0*gc27*((-mdq + muq)*sinbeta*vh + 4*cosbeta*fpi*(gpGG*(mdq + msq + muq) + (gpdd + gpss + gpuu)*vh)), 0}, {0, 0}, {(3*I)*gc27*sinbeta*vh, 0}, {(-3*I)*gc27*sinbeta*vh, 0}, {(-3*I)*gc27*sinbeta*vh, 0}, {(3*I)*gc27*sinbeta*vh, 0}, {0, 0}},

C[ S[3] , S[4] , S[1] , -S[2] ] == {{I*b0*gc28*(cosbeta*(mdq - muq)*vh + 4*fpi*sinbeta*(gpGG*(mdq + msq + muq) + (gpdd + gpss + gpuu)*vh)), 0}, {0, 0}, {(3*I)*cosbeta*gc28*vh, 0}, {(-3*I)*cosbeta*gc28*vh, 0}, {(-3*I)*cosbeta*gc28*vh, 0}, {(3*I)*cosbeta*gc28*vh, 0}, {0, 0}},

C[ S[5] , S[5] , S[2] , -S[2] ] == {{I*gc29, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[2] , -S[2] ] == {{I*b0*gc30*(2*mdq + msq + muq), 0}, {I*gc30, 0}, {(-2*I)*gc30, 0}, {I*gc30, 0}, {I*gc30, 0}, {(-2*I)*gc30, 0}, {I*gc30, 0}},

C[ S[4] , -S[4] , S[2] , -S[2] ] == {{I*b0*gc31*(mdq + msq + 2*muq), 0}, {I*gc31, 0}, {(-2*I)*gc31, 0}, {I*gc31, 0}, {I*gc31, 0}, {(-2*I)*gc31, 0}, {I*gc31, 0}},

C[ S[5] , S[7] , S[2] , -S[2] ] == {{I*gc32, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[1] , S[2] , -S[2] ] == {{I*gc33, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[7] , S[1] , S[2] , -S[2] ] == {{(-I)*b0*gc34*(-(cosbeta*(mdq + muq)*sinbeta*vh) + cosbeta^2*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh) + fpi*sinbeta^2*(gpGG*(-mdq + muq) + (-gpdd + gpuu)*vh)), 0}, {(2*I)*cosbeta*gc34*sinbeta*vh, 0}, {(-I)*cosbeta*gc34*sinbeta*vh, 0}, {(-I)*cosbeta*gc34*sinbeta*vh, 0}, {(-I)*cosbeta*gc34*sinbeta*vh, 0}, {(-I)*cosbeta*gc34*sinbeta*vh, 0}, {(2*I)*cosbeta*gc34*sinbeta*vh, 0}},

C[ S[1] , S[1] , S[2] , -S[2] ] == {{I*b0*gc35*(cosbeta*(mdq + muq)*vh + 2*fpi*sinbeta*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh)), 0}, {(2*I)*cosbeta*gc35*vh, 0}, {(-I)*cosbeta*gc35*vh, 0}, {(-I)*cosbeta*gc35*vh, 0}, {(-I)*cosbeta*gc35*vh, 0}, {(-I)*cosbeta*gc35*vh, 0}, {(2*I)*cosbeta*gc35*vh, 0}},

C[ S[2] , S[2] , -S[2] , -S[2] ] == {{(-2*I)*b0*gc36*(mdq + muq), 0}, {(2*I)*gc36, 0}, {(-I)*gc36, 0}, {(-I)*gc36, 0}, {(-I)*gc36, 0}, {(-I)*gc36, 0}, {(2*I)*gc36, 0}},

C[ S[5] , S[5] , S[7] , S[7] ] == {{I*gc37, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[7] , S[7] ] == {{(2*I)*b0*gc38*((3*mdq + msq)*sinbeta*vh - 4*cosbeta*fpi*(gpGG*(2*mdq + msq) + (2*gpdd + gpss)*vh)), 0}, {(2*I)*gc38*sinbeta*vh, 0}, {(-I)*gc38*sinbeta*vh, 0}, {(-I)*gc38*sinbeta*vh, 0}, {(-I)*gc38*sinbeta*vh, 0}, {(-I)*gc38*sinbeta*vh, 0}, {(2*I)*gc38*sinbeta*vh, 0}},

C[ S[4] , -S[4] , S[7] , S[7] ] == {{(2*I)*b0*gc39*((msq + 3*muq)*sinbeta*vh + 4*cosbeta*fpi*(gpGG*(msq + 2*muq) + (gpss + 2*gpuu)*vh)), 0}, {(2*I)*gc39*sinbeta*vh, 0}, {(-I)*gc39*sinbeta*vh, 0}, {(-I)*gc39*sinbeta*vh, 0}, {(-I)*gc39*sinbeta*vh, 0}, {(-I)*gc39*sinbeta*vh, 0}, {(2*I)*gc39*sinbeta*vh, 0}},

C[ S[5] , S[7] , S[7] , S[1] ] == {{I*gc40, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[7] , S[7] , S[1] , S[1] ] == {{I*gc41, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[7] , S[7] , S[2] , -S[2] ] == {{I*b0*gc42*((mdq + muq)*sinbeta*vh - 2*cosbeta*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh)), 0}, {(2*I)*gc42*sinbeta*vh, 0}, {(-I)*gc42*sinbeta*vh, 0}, {(-I)*gc42*sinbeta*vh, 0}, {(-I)*gc42*sinbeta*vh, 0}, {(-I)*gc42*sinbeta*vh, 0}, {(2*I)*gc42*sinbeta*vh, 0}},

C[ S[5] , S[7] , S[7] , S[7] ] == {{I*gc43, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[7] , S[7] , S[7] , S[1] ] == {{I*gc44, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[7] , S[7] , S[7] , S[7] ] == {{I*gc45, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[4] , -S[4] , V[1] , V[1] ] == {{(2*I)*qe^2, 0}},

C[ S[2] , -S[2] , V[1] , V[1] ] == {{(2*I)*qe^2, 0}},

C[ S[5] , -S[3] , -S[4] , S[2] , V[1] ] == {{(3*I)*gc48, 0}, {(-3*I)*gc48, 0}, {(-I)*gc48, 0}, {I*gc48, 0}},

C[ S[5] , S[3] , S[4] , -S[2] , V[1] ] == {{(3*I)*gc49, 0}, {(-3*I)*gc49, 0}, {(-I)*gc49, 0}, {I*gc49, 0}},

C[ S[3] , -S[3] , S[4] , -S[4] , V[1] ] == {{(3*I)*gc50, 0}, {(-3*I)*gc50, 0}, {(-I)*gc50, 0}, {I*gc50, 0}},

C[ S[3] , S[4] , S[1] , -S[2] , V[1] ] == {{I*gc51, 0}, {(-I)*gc51, 0}, {(-I)*gc51, 0}, {I*gc51, 0}},

C[ S[3] , -S[3] , S[2] , -S[2] , V[1] ] == {{(-3*I)*gc52, 0}, {(3*I)*gc52, 0}, {(-I)*gc52, 0}, {I*gc52, 0}},

C[ S[3] , S[4] , S[7] , -S[2] , V[1] ] == {{I*gc53, 0}, {(-I)*gc53, 0}, {(-I)*gc53, 0}, {I*gc53, 0}},

C[ -S[3] , -S[4] , S[1] , S[2] , V[1] ] == {{I*gc54, 0}, {(-I)*gc54, 0}, {(-I)*gc54, 0}, {I*gc54, 0}},

C[ -S[3] , -S[4] , S[7] , S[2] , V[1] ] == {{I*gc55, 0}, {(-I)*gc55, 0}, {(-I)*gc55, 0}, {I*gc55, 0}},

C[ S[4] , -S[4] , V[1] ] == {{(-I)*gc56, 0}, {I*gc56, 0}},

C[ S[5] , S[5] , S[4] , -S[4] , V[1] ] == {{0, 0}, {0, 0}, {(-I)*gc57, 0}, {I*gc57, 0}},

C[ S[4] , S[4] , -S[4] , -S[4] , V[1] ] == {{(-I)*gc58, 0}, {(-I)*gc58, 0}, {I*gc58, 0}, {I*gc58, 0}},

C[ S[5] , S[4] , -S[4] , S[1] , V[1] ] == {{0, 0}, {(-I)*gc59, 0}, {I*gc59, 0}, {0, 0}},

C[ S[4] , -S[4] , S[1] , S[1] , V[1] ] == {{(-I)*gc60, 0}, {I*gc60, 0}, {0, 0}, {0, 0}},

C[ S[4] , -S[4] , S[2] , -S[2] , V[1] ] == {{(-I)*gc61, 0}, {I*gc61, 0}, {(-I)*gc61, 0}, {I*gc61, 0}},

C[ S[5] , S[4] , -S[4] , S[7] , V[1] ] == {{0, 0}, {(-I)*gc62, 0}, {I*gc62, 0}, {0, 0}},

C[ S[4] , -S[4] , S[7] , S[1] , V[1] ] == {{(-I)*gc63, 0}, {I*gc63, 0}, {0, 0}, {0, 0}},

C[ S[4] , -S[4] , S[7] , S[7] , V[1] ] == {{(-I)*gc64, 0}, {I*gc64, 0}, {0, 0}, {0, 0}},

C[ S[2] , -S[2] , V[1] ] == {{(-I)*gc65, 0}, {I*gc65, 0}},

C[ S[1] , S[1] , S[2] , -S[2] , V[1] ] == {{0, 0}, {0, 0}, {(-I)*gc66, 0}, {I*gc66, 0}},

C[ S[2] , S[2] , -S[2] , -S[2] , V[1] ] == {{(-I)*gc67, 0}, {(-I)*gc67, 0}, {I*gc67, 0}, {I*gc67, 0}},

C[ S[7] , S[1] , S[2] , -S[2] , V[1] ] == {{0, 0}, {0, 0}, {(-I)*gc68, 0}, {I*gc68, 0}},

C[ S[7] , S[7] , S[2] , -S[2] , V[1] ] == {{0, 0}, {0, 0}, {(-I)*gc69, 0}, {I*gc69, 0}},

C[ S[5] , V[1] , V[1] ] == {{(-I)*gc70, 0}, {I*gc70, 0}, {0, 0}},

C[ S[7] , V[1] , V[1] ] == {{(-I)*gc71*sinbeta*vh, 0}, {I*gc71*sinbeta*vh, 0}, {(2*I)*cosbeta*fpi*gc71*gpFF, 0}},

C[ S[1] , V[1] , V[1] ] == {{I*cosbeta*gc72*vh, 0}, {(-I)*cosbeta*gc72*vh, 0}, {(2*I)*fpi*gc72*gpFF*sinbeta, 0}},

C[ -F[2] , F[2] , V[1] ] == {{I*gc73, 0}, {I*gc73, 0}},

C[ -F[2] , F[2] , S[7] ] == {{gc74L, 0}, {gc74R, 0}},

C[ -F[2] , F[2] , S[1] ] == {{gc75L, 0}, {gc75R, 0}},

C[ -F[3] , F[3] , S[7] ] == {{gc76L, 0}, {gc76R, 0}},

C[ -F[3] , F[3] , S[1] ] == {{gc77L, 0}, {gc77R, 0}}

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
     gc7 -> -(b0*((-mdq + muq)*sinbeta*vh + cosbeta*fpi*(gpGG*(mdq - 8*msq + muq) + (gpdd - 8*gpss + gpuu)*vh)))/(3*Sqrt[3]*fpi^2*vh),
     gc8 -> -1/(12*Sqrt[3]*fpi^2*vh),
     gc9 -> 1/(12*Sqrt[3]*fpi^2*vh),
     gc10 -> -(b0*(cosbeta*(mdq - muq)*vh + fpi*sinbeta*(gpGG*(mdq - 8*msq + muq) + (gpdd - 8*gpss + gpuu)*vh)))/(3*Sqrt[3]*fpi^2*vh),
     gc11 -> 1/(12*Sqrt[3]*fpi^2*vh),
     gc12 -> -1/(12*Sqrt[3]*fpi^2*vh),
     gc13 -> (b0*(-(cosbeta*(mdq + muq)*sinbeta*vh) + cosbeta^2*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh) + fpi*sinbeta^2*(gpGG*(-mdq + muq) + (-gpdd + gpuu)*vh)))/(3*fpi^2*vh),
     gc14 -> -1/(12*fpi^2*vh),
     gc15 -> -1/(12*fpi^2*vh),
     gc16 -> (b0*cosbeta*(cosbeta*(mdq + muq)*vh + 2*fpi*sinbeta*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh)))/(3*fpi^2*vh),
     gc17 -> cosbeta/(12*fpi^2*vh),
     gc18 -> cosbeta/(12*fpi^2*vh),
     gc19 -> -((b0*cosbeta*(cosbeta*(-mdq + muq)*sinbeta*vh + cosbeta^2*fpi*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh) - 2*fpi*sinbeta^2*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh)))/(Sqrt[3]*fpi^2*vh)),
     gc20 -> -((b0*cosbeta^2*(cosbeta*(mdq - muq)*vh + 3*fpi*sinbeta*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh)))/(Sqrt[3]*fpi^2*vh)),
     gc21 -> (b0*cosbeta^2*(-(cosbeta*(mdq + muq)*sinbeta*vh) + cosbeta^2*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh) + 3*fpi*sinbeta^2*(gpGG*(-mdq + muq) + (-gpdd + gpuu)*vh)))/(fpi^2*vh),
     gc22 -> (b0*cosbeta^3*(cosbeta*(mdq + muq)*vh + 4*fpi*sinbeta*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh)))/(fpi^2*vh),
     gc23 -> -1/(6*Sqrt[6]*fpi^2),
     gc24 -> 1/(6*Sqrt[2]*fpi^2*vh),
     gc25 -> -1/(6*Sqrt[2]*fpi^2*vh),
     gc26 -> -1/(6*Sqrt[6]*fpi^2),
     gc27 -> 1/(6*Sqrt[2]*fpi^2*vh),
     gc28 -> -1/(6*Sqrt[2]*fpi^2*vh),
     gc29 -> (b0*(mdq + muq))/(3*fpi^2),
     gc30 -> 1/(6*fpi^2),
     gc31 -> 1/(6*fpi^2),
     gc32 -> -(b0*((-mdq + muq)*sinbeta*vh + 3*cosbeta*fpi*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh)))/(3*Sqrt[3]*fpi^2*vh),
     gc33 -> -(b0*(cosbeta*(mdq - muq)*vh + 3*fpi*sinbeta*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh)))/(3*Sqrt[3]*fpi^2*vh),
     gc34 -> -1/(3*fpi^2*vh),
     gc35 -> cosbeta/(3*fpi^2*vh),
     gc36 -> -1/(3*fpi^2),
     gc37 -> (b0*sinbeta*((mdq + muq)*sinbeta*vh - 2*cosbeta*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh)))/(3*fpi^2*vh),
     gc38 -> sinbeta/(12*fpi^2*vh),
     gc39 -> sinbeta/(12*fpi^2*vh),
     gc40 -> -((b0*sinbeta*(cosbeta*(mdq - muq)*sinbeta*vh - 2*cosbeta^2*fpi*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh) + fpi*sinbeta^2*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh)))/(Sqrt[3]*fpi^2*vh)),
     gc41 -> -((b0*cosbeta*sinbeta*(-(cosbeta*(mdq + muq)*sinbeta*vh) + 2*cosbeta^2*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh) + 2*fpi*sinbeta^2*(gpGG*(-mdq + muq) + (-gpdd + gpuu)*vh)))/(fpi^2*vh)),
     gc42 -> sinbeta/(3*fpi^2*vh),
     gc43 -> -((b0*sinbeta^2*((-mdq + muq)*sinbeta*vh + 3*cosbeta*fpi*(gpGG*(mdq + muq) + (gpdd + gpuu)*vh)))/(Sqrt[3]*fpi^2*vh)),
     gc44 -> -((b0*sinbeta^2*(cosbeta*(mdq + muq)*sinbeta*vh - 3*cosbeta^2*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh) + fpi*sinbeta^2*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh)))/(fpi^2*vh)),
     gc45 -> (b0*sinbeta^3*((mdq + muq)*sinbeta*vh - 4*cosbeta*fpi*(gpGG*(mdq - muq) + (gpdd - gpuu)*vh)))/(fpi^2*vh),
     gc48 -> -qe/(2*Sqrt[6]*fpi^2),
     gc49 -> qe/(2*Sqrt[6]*fpi^2),
     gc50 -> -qe/(6*fpi^2),
     gc51 -> (cosbeta*qe)/(2*Sqrt[2]*fpi^2),
     gc52 -> -qe/(6*fpi^2),
     gc53 -> -(qe*sinbeta)/(2*Sqrt[2]*fpi^2),
     gc54 -> -(cosbeta*qe)/(2*Sqrt[2]*fpi^2),
     gc55 -> (qe*sinbeta)/(2*Sqrt[2]*fpi^2),
     gc56 -> qe,
     gc57 -> -qe/(2*fpi^2),
     gc58 -> (-4*qe)/(3*fpi^2),
     gc59 -> -(cosbeta*qe)/(2*Sqrt[3]*fpi^2),
     gc60 -> -(cosbeta^2*qe)/(6*fpi^2),
     gc61 -> (-2*qe)/(3*fpi^2),
     gc62 -> (qe*sinbeta)/(2*Sqrt[3]*fpi^2),
     gc63 -> (cosbeta*qe*sinbeta)/(6*fpi^2),
     gc64 -> -(qe*sinbeta^2)/(6*fpi^2),
     gc65 -> qe,
     gc66 -> (-2*cosbeta^2*qe)/(3*fpi^2),
     gc67 -> (-4*qe)/(3*fpi^2),
     gc68 -> (2*cosbeta*qe*sinbeta)/(3*fpi^2),
     gc69 -> (-2*qe*sinbeta^2)/(3*fpi^2),
     gc70 -> -alphaEM/(2*Sqrt[3]*fpi*Pi),
     gc71 -> alphaEM/(2*fpi*Pi*vh),
     gc72 -> alphaEM/(2*fpi*Pi*vh),
     gc73 -> -qe,
     gc74L -> cosbeta*gpll,
     gc74R -> -(cosbeta*gpll),
     gc75L -> gpll*sinbeta,
     gc75R -> -(gpll*sinbeta),
     gc76L -> cosbeta*gpxx,
     gc76R -> -(cosbeta*gpxx),
     gc77L -> gpxx*sinbeta,
     gc77R -> -(gpxx*sinbeta)};
