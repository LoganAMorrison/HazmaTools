(* Patched for use with FeynCalc *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"Hazma",
  Authors -> {"Adam Coogan", "Logan Morrison"},
  Date -> "11/20/2019",
  Institutions -> {"University of Amsterdam, GRAPPA", "University of California, Santa Cruz", "Santa Cruz Institute for Particle Physics"},
  Emails -> {"a.m.coogan@uva.nl", "loanmorr@ucsc.edu"},
  URLs -> {"https://github.com/LoganAMorrison/HazmaTools", "https://github.com/LoganAMorrison/Hazma"},
  References -> {"https://arxiv.org/pdf/1907.11846.pdf"},
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
    PropagatorLabel -> "S",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> ms,
    Indices -> {} },

S[8] == {
    SelfConjugate -> True,
    PropagatorLabel -> "P",
    PropagatorType -> ScalarDash,
    PropagatorArrow -> None,
    Mass -> mp,
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
    PropagatorLabel -> "\\rho",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> mrho,
    Indices -> {} },

V[3] == {
    SelfConjugate -> True,
    PropagatorLabel -> "\\omega",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> momega,
    Indices -> {} },

V[4] == {
    SelfConjugate -> True,
    PropagatorLabel -> "V",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> mv,
    Indices -> {} },

V[5] == {
    SelfConjugate -> True,
    PropagatorLabel -> "A",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> ma,
    Indices -> {} },

F[1] == {
    SelfConjugate -> False,
    QuantumNumbers -> {LeptonNumber},
    PropagatorLabel -> ComposedChar["\\nu", "\\ell", Null],
    PropagatorType -> Straight,
    PropagatorArrow -> Forward,
    Mass -> 0,
    Indices -> {} },

F[2] == {
    SelfConjugate -> False,
    QuantumNumbers -> {-Q, LeptonNumber},
    PropagatorLabel -> "\\ell",
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

FAGaugeXi[ V[1] ] = FAGaugeXi[A];

mpi0[ ___ ] := mpi0;
mpi[ ___ ] := mpi;
mk0[ ___ ] := mk0;
mk[ ___ ] := mk;
meta[ ___ ] := meta;
ms[ ___ ] := ms;
mp[ ___ ] := mp;
mrho[ ___ ] := mrho;
momega[ ___ ] := momega;
mv[ ___ ] := mv;
ma[ ___ ] := ma;
ml[ ___ ] := ml;
mx[ ___ ] := mx;




(*      Couplings (calculated by FeynRules)      *)

M$CouplingMatrices = {

C[ S[5] , S[5] , S[5] , S[5] ] == {{((I/9)*b0*(mdq + 16*msq + muq))/fpi^2, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[5] , S[3] , -S[3] ] == {{((I/6)*b0*(mdq + 3*msq))/fpi^2, 0}, {(I/2)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(I/2)/fpi^2, 0}},

C[ S[3] , S[3] , -S[3] , -S[3] ] == {{(((2*I)/3)*b0*(mdq + msq))/fpi^2, 0}, {((-2*I)/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {((-2*I)/3)/fpi^2, 0}},

C[ S[3] , -S[3] , S[4] , -S[4] ] == {{((I/6)*b0*(mdq + 2*msq + muq))/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(I/6)/fpi^2, 0}},

C[ S[5] , S[5] , S[4] , -S[4] ] == {{((I/6)*b0*(3*msq + muq))/fpi^2, 0}, {(I/2)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(-I/4)/fpi^2, 0}, {(I/2)/fpi^2, 0}},

C[ S[4] , S[4] , -S[4] , -S[4] ] == {{(((2*I)/3)*b0*(msq + muq))/fpi^2, 0}, {((-2*I)/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {((-2*I)/3)/fpi^2, 0}},

C[ S[5] , S[5] , S[5] , S[1] ] == {{((-I/3)*b0*(mdq - muq))/(Sqrt[3]*fpi^2), 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[3] , -S[3] , S[1] ] == {{((-I/6)*b0*(mdq - msq))/(Sqrt[3]*fpi^2), 0}, {(I/4)/(Sqrt[3]*fpi^2), 0}, {(I/4)/(Sqrt[3]*fpi^2), 0}, {(-I/2)/(Sqrt[3]*fpi^2), 0}, {(-I/2)/(Sqrt[3]*fpi^2), 0}, {(I/4)/(Sqrt[3]*fpi^2), 0}, {(I/4)/(Sqrt[3]*fpi^2), 0}},

C[ S[5] , S[4] , -S[4] , S[1] ] == {{((-I/6)*b0*(msq - muq))/(Sqrt[3]*fpi^2), 0}, {(-I/4)/(Sqrt[3]*fpi^2), 0}, {(-I/4)/(Sqrt[3]*fpi^2), 0}, {(I/2)/(Sqrt[3]*fpi^2), 0}, {(I/2)/(Sqrt[3]*fpi^2), 0}, {(-I/4)/(Sqrt[3]*fpi^2), 0}, {(-I/4)/(Sqrt[3]*fpi^2), 0}},

C[ S[5] , S[5] , S[1] , S[1] ] == {{((I/3)*b0*(mdq + muq))/fpi^2, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[1] , S[1] ] == {{((I/6)*b0*(3*mdq + msq))/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(I/6)/fpi^2, 0}},

C[ S[4] , -S[4] , S[1] , S[1] ] == {{((I/6)*b0*(msq + 3*muq))/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(-I/12)/fpi^2, 0}, {(I/6)/fpi^2, 0}},

C[ S[5] , S[1] , S[1] , S[1] ] == {{((-I)*b0*(mdq - muq))/(Sqrt[3]*fpi^2), 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[1] , S[1] , S[1] , S[1] ] == {{(I*b0*(mdq + muq))/fpi^2, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , -S[3] , -S[4] , S[2] ] == {{((I/6)*b0*(mdq - 2*msq + muq))/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}, {I/(Sqrt[6]*fpi^2), 0}, {I/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}},

C[ -S[3] , -S[4] , S[1] , S[2] ] == {{((-I/6)*b0*(mdq - muq))/(Sqrt[2]*fpi^2), 0}, {0, 0}, {(-I/2)/(Sqrt[2]*fpi^2), 0}, {(I/2)/(Sqrt[2]*fpi^2), 0}, {(I/2)/(Sqrt[2]*fpi^2), 0}, {(-I/2)/(Sqrt[2]*fpi^2), 0}, {0, 0}},

C[ S[5] , S[3] , S[4] , -S[2] ] == {{((I/6)*b0*(mdq - 2*msq + muq))/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}, {I/(Sqrt[6]*fpi^2), 0}, {I/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}, {(-I/2)/(Sqrt[6]*fpi^2), 0}},

C[ S[3] , S[4] , S[1] , -S[2] ] == {{((-I/6)*b0*(mdq - muq))/(Sqrt[2]*fpi^2), 0}, {0, 0}, {(-I/2)/(Sqrt[2]*fpi^2), 0}, {(I/2)/(Sqrt[2]*fpi^2), 0}, {(I/2)/(Sqrt[2]*fpi^2), 0}, {(-I/2)/(Sqrt[2]*fpi^2), 0}, {0, 0}},

C[ S[5] , S[5] , S[2] , -S[2] ] == {{((I/3)*b0*(mdq + muq))/fpi^2, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[2] , -S[2] ] == {{((I/6)*b0*(2*mdq + msq + muq))/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(I/6)/fpi^2, 0}},

C[ S[4] , -S[4] , S[2] , -S[2] ] == {{((I/6)*b0*(mdq + msq + 2*muq))/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(I/6)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(I/6)/fpi^2, 0}},

C[ S[5] , S[1] , S[2] , -S[2] ] == {{((-I/3)*b0*(mdq - muq))/(Sqrt[3]*fpi^2), 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[1] , S[1] , S[2] , -S[2] ] == {{((I/3)*b0*(mdq + muq))/fpi^2, 0}, {((2*I)/3)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {(-I/3)/fpi^2, 0}, {((2*I)/3)/fpi^2, 0}},

C[ S[2] , S[2] , -S[2] , -S[2] ] == {{(((2*I)/3)*b0*(mdq + muq))/fpi^2, 0}, {((-2*I)/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {(I/3)/fpi^2, 0}, {((-2*I)/3)/fpi^2, 0}},

C[ S[4] , -S[4] , V[1] , V[1] ] == {{(2*I)*qe^2, 0}},

C[ S[2] , -S[2] , V[1] , V[1] ] == {{(2*I)*qe^2, 0}},

C[ S[5] , -S[3] , -S[4] , S[2] , V[1] ] == {{((-I/2)*Sqrt[3/2]*qe)/fpi^2, 0}, {((I/2)*Sqrt[3/2]*qe)/fpi^2, 0}, {((I/2)*qe)/(Sqrt[6]*fpi^2), 0}, {((-I/2)*qe)/(Sqrt[6]*fpi^2), 0}},

C[ S[5] , S[3] , S[4] , -S[2] , V[1] ] == {{((I/2)*Sqrt[3/2]*qe)/fpi^2, 0}, {((-I/2)*Sqrt[3/2]*qe)/fpi^2, 0}, {((-I/2)*qe)/(Sqrt[6]*fpi^2), 0}, {((I/2)*qe)/(Sqrt[6]*fpi^2), 0}},

C[ S[3] , -S[3] , S[4] , -S[4] , V[1] ] == {{((-I/2)*qe)/fpi^2, 0}, {((I/2)*qe)/fpi^2, 0}, {((I/6)*qe)/fpi^2, 0}, {((-I/6)*qe)/fpi^2, 0}},

C[ S[3] , S[4] , S[1] , -S[2] , V[1] ] == {{((I/2)*qe)/(Sqrt[2]*fpi^2), 0}, {((-I/2)*qe)/(Sqrt[2]*fpi^2), 0}, {((-I/2)*qe)/(Sqrt[2]*fpi^2), 0}, {((I/2)*qe)/(Sqrt[2]*fpi^2), 0}},

C[ S[3] , -S[3] , S[2] , -S[2] , V[1] ] == {{((I/2)*qe)/fpi^2, 0}, {((-I/2)*qe)/fpi^2, 0}, {((I/6)*qe)/fpi^2, 0}, {((-I/6)*qe)/fpi^2, 0}},

C[ -S[3] , -S[4] , S[1] , S[2] , V[1] ] == {{((-I/2)*qe)/(Sqrt[2]*fpi^2), 0}, {((I/2)*qe)/(Sqrt[2]*fpi^2), 0}, {((I/2)*qe)/(Sqrt[2]*fpi^2), 0}, {((-I/2)*qe)/(Sqrt[2]*fpi^2), 0}},

C[ S[4] , -S[4] , V[1] ] == {{(-I)*qe, 0}, {I*qe, 0}},

C[ S[5] , S[5] , S[4] , -S[4] , V[1] ] == {{0, 0}, {0, 0}, {((I/2)*qe)/fpi^2, 0}, {((-I/2)*qe)/fpi^2, 0}},

C[ S[4] , S[4] , -S[4] , -S[4] , V[1] ] == {{(((4*I)/3)*qe)/fpi^2, 0}, {(((4*I)/3)*qe)/fpi^2, 0}, {(((-4*I)/3)*qe)/fpi^2, 0}, {(((-4*I)/3)*qe)/fpi^2, 0}},

C[ S[5] , S[4] , -S[4] , S[1] , V[1] ] == {{0, 0}, {((I/2)*qe)/(Sqrt[3]*fpi^2), 0}, {((-I/2)*qe)/(Sqrt[3]*fpi^2), 0}, {0, 0}},

C[ S[4] , -S[4] , S[1] , S[1] , V[1] ] == {{((I/6)*qe)/fpi^2, 0}, {((-I/6)*qe)/fpi^2, 0}, {0, 0}, {0, 0}},

C[ S[4] , -S[4] , S[2] , -S[2] , V[1] ] == {{(((2*I)/3)*qe)/fpi^2, 0}, {(((-2*I)/3)*qe)/fpi^2, 0}, {(((2*I)/3)*qe)/fpi^2, 0}, {(((-2*I)/3)*qe)/fpi^2, 0}},

C[ S[2] , -S[2] , V[1] ] == {{(-I)*qe, 0}, {I*qe, 0}},

C[ S[1] , S[1] , S[2] , -S[2] , V[1] ] == {{0, 0}, {0, 0}, {(((2*I)/3)*qe)/fpi^2, 0}, {(((-2*I)/3)*qe)/fpi^2, 0}},

C[ S[2] , S[2] , -S[2] , -S[2] , V[1] ] == {{(((4*I)/3)*qe)/fpi^2, 0}, {(((4*I)/3)*qe)/fpi^2, 0}, {(((-4*I)/3)*qe)/fpi^2, 0}, {(((-4*I)/3)*qe)/fpi^2, 0}},

C[ S[5] , S[3] , -S[3] , V[1] ] == {{0, 0}, {((-I/4)*Sqrt[3]*qe)/(fpi^3*Pi^2), 0}},

C[ S[5] , S[4] , -S[4] , V[1] , V[1] ] == {{((-I/8)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , S[4] , -S[2] , V[1] , V[1] ] == {{((-I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -S[3] , -S[4] , S[2] , V[1] , V[1] ] == {{((-I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(Sqrt[2]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[4] , -S[4] , V[1] ] == {{0, 0}, {((-I/4)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}},

C[ S[3] , -S[3] , S[1] , V[1] ] == {{0, 0}, {((I/4)*qe)/(fpi^3*Pi^2), 0}},

C[ S[4] , -S[4] , S[1] , V[1] , V[1] ] == {{((I/24)*qe^2)/(fpi^3*Pi^2), 0}, {((-I/24)*qe^2)/(fpi^3*Pi^2), 0}, {((I/24)*qe^2)/(fpi^3*Pi^2), 0}, {((-I/24)*qe^2)/(fpi^3*Pi^2), 0}, {(((-5*I)/24)*qe^2)/(fpi^3*Pi^2), 0}, {(((5*I)/24)*qe^2)/(fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[4] , -S[4] , S[1] , V[1] ] == {{0, 0}, {((-I/4)*qe)/(fpi^3*Pi^2), 0}},

C[ S[5] , S[2] , -S[2] , V[1] , V[1] ] == {{((-I/4)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/4)*qe^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[1] , S[2] , -S[2] , V[1] , V[1] ] == {{((-I/6)*qe^2)/(fpi^3*Pi^2), 0}, {((I/6)*qe^2)/(fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(fpi^3*Pi^2), 0}, {((-I/12)*qe^2)/(fpi^3*Pi^2), 0}, {((I/12)*qe^2)/(fpi^3*Pi^2), 0}, {((-I/12)*qe^2)/(fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[2] , -S[2] , V[1] ] == {{0, 0}, {((-I/4)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}},

C[ S[1] , S[2] , -S[2] , V[1] ] == {{0, 0}, {((-I/4)*qe)/(fpi^3*Pi^2), 0}},

C[ S[5] , V[1] , V[1] ] == {{((I/8)*qe^2)/(Sqrt[3]*fpi*Pi^2), 0}, {((-I/8)*qe^2)/(Sqrt[3]*fpi*Pi^2), 0}},

C[ S[1] , V[1] , V[1] ] == {{((I/8)*qe^2)/(fpi*Pi^2), 0}, {((-I/8)*qe^2)/(fpi*Pi^2), 0}},

C[ S[5] , V[1] , V[4] ] == {{((-I/8)*(gvdd - 2*(gvss + gvuu))*qe)/(Sqrt[3]*fpi*Pi^2), 0}, {((I/8)*(gvdd - 2*(gvss + gvuu))*qe)/(Sqrt[3]*fpi*Pi^2), 0}},

C[ S[5] , S[3] , -S[3] , V[1] , V[4] ] == {{((I/4)*(gvdd - gvss)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/4)*(gvdd - gvss)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*(gvdd - gvss)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*(gvdd - gvss)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*(gvdd - gvss)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*(gvdd - gvss)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {0, 0}},

C[ S[5] , S[4] , -S[4] , V[1] , V[4] ] == {{((-I/8)*Sqrt[3]*(gvss + gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/8)*Sqrt[3]*(gvss + gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/8)*(gvss + 2*gvuu)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*(gvss + 2*gvuu)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*(gvss + 2*gvuu)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*(gvss + 2*gvuu)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(2*gvss + gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(2*gvss + gvuu)*qe)/(fpi^3*Pi^2), 0}, {0, 0}},

C[ S[5] , S[2] , -S[2] , V[1] , V[4] ] == {{((I/4)*(gvdd - gvuu)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/4)*(gvdd - gvuu)*qe)/(Sqrt[3]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[1] , V[1] , V[4] ] == {{((I/24)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {((-I/24)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {((I/24)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {((-I/24)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {((-I/12)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {((I/12)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {0, 0}, {((-I/8)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}, {((-I/8)*(gvdd - gvss)*qe)/(fpi^3*Pi^2), 0}},

C[ S[3] , S[4] , -S[2] , V[1] , V[4] ] == {{((I/24)*(5*gvdd + 5*gvss + 2*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/24)*(5*gvdd + 5*gvss + 2*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd + gvss + gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd + gvss + gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {0, 0}},

C[ -S[3] , -S[4] , S[2] , V[1] , V[4] ] == {{((I/24)*(5*gvdd + 5*gvss + 2*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/24)*(5*gvdd + 5*gvss + 2*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/24)*(gvdd + gvss + 4*gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd + gvss + gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd + gvss + gvuu)*qe)/(Sqrt[2]*fpi^3*Pi^2), 0}, {0, 0}},

C[ S[4] , -S[4] , S[1] , V[1] , V[4] ] == {{((I/24)*(gvss + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/24)*(gvss + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/24)*(gvss + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/24)*(gvss + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/24)*(gvss - 7*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/24)*(gvss - 7*gvuu)*qe)/(fpi^3*Pi^2), 0}, {0, 0}, {((I/8)*(2*gvss + gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/8)*(2*gvss + gvuu)*qe)/(fpi^3*Pi^2), 0}},

C[ S[1] , V[1] , V[4] ] == {{((I/8)*(gvdd + 2*gvuu)*qe)/(fpi*Pi^2), 0}, {((-I/8)*(gvdd + 2*gvuu)*qe)/(fpi*Pi^2), 0}},

C[ S[1] , S[2] , -S[2] , V[1] , V[4] ] == {{((-I/6)*(gvdd + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/6)*(gvdd + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/12)*(gvdd + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/12)*(gvdd + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((I/12)*(gvdd + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/12)*(gvdd + 2*gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/4)*(2*gvdd + gvuu)*qe)/(fpi^3*Pi^2), 0}, {((-I/4)*(2*gvdd + gvuu)*qe)/(fpi^3*Pi^2), 0}, {0, 0}},

C[ -F[2] , F[2] , V[1] ] == {{(-I)*qe, 0}, {(-I)*qe, 0}},

C[ -F[2] , F[2] , V[4] ] == {{I*gvll, 0}, {I*gvll, 0}},

C[ -F[3] , F[3] , V[4] ] == {{I*gvxx, 0}, {I*gvxx, 0}},

C[ S[4] , -S[4] , V[1] , V[4] ] == {{(-2*I)*(gvss - gvuu)*qe, 0}},

C[ S[2] , -S[2] , V[1] , V[4] ] == {{(-2*I)*(gvdd - gvuu)*qe, 0}},

C[ S[5] , -S[3] , -S[4] , S[2] , V[4] ] == {{((I/2)*Sqrt[3/2]*(gvdd - gvuu))/fpi^2, 0}, {((-I/2)*(gvdd + 2*gvss - 3*gvuu))/(Sqrt[6]*fpi^2), 0}, {((-I/2)*(3*gvdd - 2*gvss - gvuu))/(Sqrt[6]*fpi^2), 0}, {((I/2)*(gvdd - gvuu))/(Sqrt[6]*fpi^2), 0}},

C[ S[5] , S[3] , S[4] , -S[2] , V[4] ] == {{((-I/2)*Sqrt[3/2]*(gvdd - gvuu))/fpi^2, 0}, {((I/2)*(gvdd + 2*gvss - 3*gvuu))/(Sqrt[6]*fpi^2), 0}, {((I/2)*(3*gvdd - 2*gvss - gvuu))/(Sqrt[6]*fpi^2), 0}, {((-I/2)*(gvdd - gvuu))/(Sqrt[6]*fpi^2), 0}},

C[ S[3] , -S[3] , V[4] ] == {{I*(gvdd - gvss), 0}, {(-I)*(gvdd - gvss), 0}},

C[ S[5] , S[5] , S[3] , -S[3] , V[4] ] == {{0, 0}, {0, 0}, {((-I/2)*(gvdd - gvss))/fpi^2, 0}, {((I/2)*(gvdd - gvss))/fpi^2, 0}},

C[ S[3] , S[3] , -S[3] , -S[3] , V[4] ] == {{(((-4*I)/3)*(gvdd - gvss))/fpi^2, 0}, {(((-4*I)/3)*(gvdd - gvss))/fpi^2, 0}, {(((4*I)/3)*(gvdd - gvss))/fpi^2, 0}, {(((4*I)/3)*(gvdd - gvss))/fpi^2, 0}},

C[ S[3] , -S[3] , S[4] , -S[4] , V[4] ] == {{((-I/6)*(gvdd - 4*gvss + 3*gvuu))/fpi^2, 0}, {((I/6)*(gvdd - 4*gvss + 3*gvuu))/fpi^2, 0}, {((I/6)*(3*gvdd - 4*gvss + gvuu))/fpi^2, 0}, {((-I/6)*(3*gvdd - 4*gvss + gvuu))/fpi^2, 0}},

C[ S[5] , S[3] , -S[3] , S[1] , V[4] ] == {{0, 0}, {((I/2)*(gvdd - gvss))/(Sqrt[3]*fpi^2), 0}, {((-I/2)*(gvdd - gvss))/(Sqrt[3]*fpi^2), 0}, {0, 0}},

C[ S[3] , -S[3] , S[1] , S[1] , V[4] ] == {{((-I/6)*(gvdd - gvss))/fpi^2, 0}, {((I/6)*(gvdd - gvss))/fpi^2, 0}, {0, 0}, {0, 0}},

C[ S[3] , S[4] , S[1] , -S[2] , V[4] ] == {{((-I/2)*(gvdd - gvuu))/(Sqrt[2]*fpi^2), 0}, {((I/2)*(gvdd - gvuu))/(Sqrt[2]*fpi^2), 0}, {((-I/2)*(gvdd - 2*gvss + gvuu))/(Sqrt[2]*fpi^2), 0}, {((I/2)*(gvdd - 2*gvss + gvuu))/(Sqrt[2]*fpi^2), 0}},

C[ S[3] , -S[3] , S[2] , -S[2] , V[4] ] == {{((-I/6)*(4*gvdd - gvss - 3*gvuu))/fpi^2, 0}, {((I/6)*(4*gvdd - gvss - 3*gvuu))/fpi^2, 0}, {((-I/6)*(4*gvdd - 3*gvss - gvuu))/fpi^2, 0}, {((I/6)*(4*gvdd - 3*gvss - gvuu))/fpi^2, 0}},

C[ -S[3] , -S[4] , S[1] , S[2] , V[4] ] == {{((I/2)*(gvdd - gvuu))/(Sqrt[2]*fpi^2), 0}, {((-I/2)*(gvdd - gvuu))/(Sqrt[2]*fpi^2), 0}, {((I/2)*(gvdd - 2*gvss + gvuu))/(Sqrt[2]*fpi^2), 0}, {((-I/2)*(gvdd - 2*gvss + gvuu))/(Sqrt[2]*fpi^2), 0}},

C[ S[4] , -S[4] , V[4] ] == {{I*(gvss - gvuu), 0}, {(-I)*(gvss - gvuu), 0}},

C[ S[5] , S[5] , S[4] , -S[4] , V[4] ] == {{0, 0}, {0, 0}, {((-I/2)*(gvss - gvuu))/fpi^2, 0}, {((I/2)*(gvss - gvuu))/fpi^2, 0}},

C[ S[4] , S[4] , -S[4] , -S[4] , V[4] ] == {{(((-4*I)/3)*(gvss - gvuu))/fpi^2, 0}, {(((-4*I)/3)*(gvss - gvuu))/fpi^2, 0}, {(((4*I)/3)*(gvss - gvuu))/fpi^2, 0}, {(((4*I)/3)*(gvss - gvuu))/fpi^2, 0}},

C[ S[5] , S[4] , -S[4] , S[1] , V[4] ] == {{0, 0}, {((-I/2)*(gvss - gvuu))/(Sqrt[3]*fpi^2), 0}, {((I/2)*(gvss - gvuu))/(Sqrt[3]*fpi^2), 0}, {0, 0}},

C[ S[4] , -S[4] , S[1] , S[1] , V[4] ] == {{((-I/6)*(gvss - gvuu))/fpi^2, 0}, {((I/6)*(gvss - gvuu))/fpi^2, 0}, {0, 0}, {0, 0}},

C[ S[4] , -S[4] , S[2] , -S[2] , V[4] ] == {{((-I/6)*(3*gvdd + gvss - 4*gvuu))/fpi^2, 0}, {((I/6)*(3*gvdd + gvss - 4*gvuu))/fpi^2, 0}, {((-I/6)*(gvdd + 3*gvss - 4*gvuu))/fpi^2, 0}, {((I/6)*(gvdd + 3*gvss - 4*gvuu))/fpi^2, 0}},

C[ S[2] , -S[2] , V[4] ] == {{I*(gvdd - gvuu), 0}, {(-I)*(gvdd - gvuu), 0}},

C[ S[1] , S[1] , S[2] , -S[2] , V[4] ] == {{0, 0}, {0, 0}, {(((-2*I)/3)*(gvdd - gvuu))/fpi^2, 0}, {(((2*I)/3)*(gvdd - gvuu))/fpi^2, 0}},

C[ S[2] , S[2] , -S[2] , -S[2] , V[4] ] == {{(((-4*I)/3)*(gvdd - gvuu))/fpi^2, 0}, {(((-4*I)/3)*(gvdd - gvuu))/fpi^2, 0}, {(((4*I)/3)*(gvdd - gvuu))/fpi^2, 0}, {(((4*I)/3)*(gvdd - gvuu))/fpi^2, 0}},

C[ S[3] , -S[3] , V[4] , V[4] ] == {{(2*I)*(gvdd - gvss)^2, 0}},

C[ S[4] , -S[4] , V[4] , V[4] ] == {{(2*I)*(gvss - gvuu)^2, 0}},

C[ S[2] , -S[2] , V[4] , V[4] ] == {{(2*I)*(gvdd - gvuu)^2, 0}},

C[ S[5] , V[4] , V[4] ] == {{((I/8)*Sqrt[3]*(gvdd^2 - 2*gvss^2 + gvuu^2))/(fpi*Pi^2), 0}, {((-I/8)*Sqrt[3]*(gvdd^2 - 2*gvss^2 + gvuu^2))/(fpi*Pi^2), 0}},

C[ S[5] , S[3] , -S[3] , V[4] , V[4] ] == {{((-I/8)*(gvdd - gvss)*(5*gvdd + 7*gvss))/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*(gvdd - gvss)*(5*gvdd + 7*gvss))/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/8)*Sqrt[3]*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {((I/8)*Sqrt[3]*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[4] , -S[4] , V[4] , V[4] ] == {{((I/8)*(gvss - gvuu)*(7*gvss + 5*gvuu))/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*(gvss - gvuu)*(7*gvss + 5*gvuu))/(Sqrt[3]*fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {((I/8)*Sqrt[3]*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {((-I/8)*Sqrt[3]*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {((I/8)*Sqrt[3]*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[2] , -S[2] , V[4] , V[4] ] == {{((-I/4)*(gvdd - gvuu)^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd - gvuu)^2)/(Sqrt[3]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , -S[3] , S[1] , V[4] , V[4] ] == {{((-I/8)*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {((I/8)*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {((-I/8)*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {((I/8)*(gvdd - gvss)*(gvdd + gvss))/(fpi^3*Pi^2), 0}, {((I/8)*(gvdd - gvss)*(3*gvdd + gvss))/(fpi^3*Pi^2), 0}, {((-I/8)*(gvdd - gvss)*(3*gvdd + gvss))/(fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[3] , S[4] , -S[2] , V[4] , V[4] ] == {{((-I/4)*(gvdd^2 + gvss^2 + gvdd*(gvss - gvuu) - gvss*gvuu - gvuu^2))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd^2 + gvss^2 + gvdd*(gvss - gvuu) - gvss*gvuu - gvuu^2))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd^2 - gvss^2 - gvss*gvuu - gvuu^2 + gvdd*(gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/4)*(gvdd^2 - gvss^2 - gvss*gvuu - gvuu^2 + gvdd*(gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/4)*(gvdd^2 - gvss^2 - gvss*gvuu + gvuu^2 + gvdd*(-gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd^2 - gvss^2 - gvss*gvuu + gvuu^2 + gvdd*(-gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -S[3] , -S[4] , S[2] , V[4] , V[4] ] == {{((-I/4)*(gvdd^2 + gvss^2 + gvdd*(gvss - gvuu) - gvss*gvuu - gvuu^2))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd^2 + gvss^2 + gvdd*(gvss - gvuu) - gvss*gvuu - gvuu^2))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd^2 - gvss^2 - gvss*gvuu - gvuu^2 + gvdd*(gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/4)*(gvdd^2 - gvss^2 - gvss*gvuu - gvuu^2 + gvdd*(gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((-I/4)*(gvdd^2 - gvss^2 - gvss*gvuu + gvuu^2 + gvdd*(-gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {((I/4)*(gvdd^2 - gvss^2 - gvss*gvuu + gvuu^2 + gvdd*(-gvss + gvuu)))/(Sqrt[2]*fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[4] , -S[4] , S[1] , V[4] , V[4] ] == {{((-I/8)*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {((I/8)*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {((-I/8)*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {((I/8)*(gvss - gvuu)*(gvss + gvuu))/(fpi^3*Pi^2), 0}, {((I/8)*(gvss - gvuu)*(gvss + 3*gvuu))/(fpi^3*Pi^2), 0}, {((-I/8)*(gvss - gvuu)*(gvss + 3*gvuu))/(fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[1] , V[4] , V[4] ] == {{(((-3*I)/8)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi*Pi^2), 0}, {(((3*I)/8)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi*Pi^2), 0}},

C[ S[1] , S[2] , -S[2] , V[4] , V[4] ] == {{((I/2)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi^3*Pi^2), 0}, {((-I/2)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi^3*Pi^2), 0}, {((-I/4)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi^3*Pi^2), 0}, {((I/4)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi^3*Pi^2), 0}, {((-I/4)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi^3*Pi^2), 0}, {((I/4)*(gvdd - gvuu)*(gvdd + gvuu))/(fpi^3*Pi^2), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ S[5] , S[3] , -S[3] , V[4] ] == {{0, 0}, {((I/4)*(4*gvdd + 5*gvss))/(Sqrt[3]*fpi^3*Pi^2), 0}},

C[ S[5] , S[4] , -S[4] , V[4] ] == {{0, 0}, {((-I/4)*(5*gvss + 4*gvuu))/(Sqrt[3]*fpi^3*Pi^2), 0}},

C[ S[3] , -S[3] , S[1] , V[4] ] == {{0, 0}, {((-I/4)*(2*gvdd + gvss))/(fpi^3*Pi^2), 0}},

C[ S[4] , -S[4] , S[1] , V[4] ] == {{0, 0}, {((-I/4)*(gvss + 2*gvuu))/(fpi^3*Pi^2), 0}},

C[ -S[3] , -S[4] , S[2] , V[4] ] == {{0, 0}, {((-I/2)*(gvdd + gvss + gvuu))/(Sqrt[2]*fpi^3*Pi^2), 0}},

C[ S[3] , S[4] , -S[2] , V[4] ] == {{0, 0}, {((I/2)*(gvdd + gvss + gvuu))/(Sqrt[2]*fpi^3*Pi^2), 0}},

C[ S[5] , S[2] , -S[2] , V[4] ] == {{0, 0}, {((I/4)*(gvdd - gvuu))/(Sqrt[3]*fpi^3*Pi^2), 0}},

C[ S[1] , S[2] , -S[2] , V[4] ] == {{0, 0}, {(((-3*I)/4)*(gvdd + gvuu))/(fpi^3*Pi^2), 0}},

C[ -F[1] , F[2] , -S[4] , S[2] , -S[2] ] == {{(-2*GF*Vus)/(3*fpi), 0}, {(4*GF*Vus)/(3*fpi), 0}, {(-2*GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[2] ] == {{2*fpi*GF*Vud, 0}, {0, 0}},

C[ -F[1] , F[2] , S[3] , -S[3] , -S[2] ] == {{(4*GF*Vud)/(3*fpi), 0}, {(-2*GF*Vud)/(3*fpi), 0}, {(-2*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[4] , -S[4] , -S[2] ] == {{(4*GF*Vud)/(3*fpi), 0}, {(-2*GF*Vud)/(3*fpi), 0}, {(-2*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[1] , -S[2] ] == {{(2*I)*GF*Vud, 0}, {(-2*I)*GF*Vud, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[1] , S[1] , -S[2] ] == {{(4*GF*Vud)/(3*fpi), 0}, {(4*GF*Vud)/(3*fpi), 0}, {(-8*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[2] , -S[2] , -S[2] ] == {{(8*GF*Vud)/(3*fpi), 0}, {(-4*GF*Vud)/(3*fpi), 0}, {(-4*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[3] , -S[2] ] == {{(-I)*Sqrt[2]*GF*Vus, 0}, {I*Sqrt[2]*GF*Vus, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[5] , S[3] , -S[2] ] == {{(Sqrt[2/3]*GF*Vus)/fpi, 0}, {(-2*Sqrt[2/3]*GF*Vus)/fpi, 0}, {(Sqrt[2/3]*GF*Vus)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[3] , S[1] , -S[2] ] == {{0, 0}, {-((Sqrt[2]*GF*Vus)/fpi), 0}, {(Sqrt[2]*GF*Vus)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[3] , S[4] , V[4] ] == {{I*Sqrt[2]*GF*(gvdd - 2*gvss + gvuu)*Vud, 0}, {0, 0}},

C[ -F[2] , F[1] , S[2] , V[4] ] == {{2*fpi*GF*(gvdd - gvuu)*Vud, 0}, {0, 0}},

C[ -F[2] , F[1] , S[1] , S[2] , V[4] ] == {{(2*I)*GF*(gvdd - gvuu)*Vud, 0}, {0, 0}},

C[ -F[2] , F[1] , S[5] , S[4] , V[4] ] == {{I*Sqrt[3]*GF*(gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , V[4] ] == {{2*fpi*GF*(gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , S[1] , V[4] ] == {{I*GF*(gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , -S[3] , S[2] , V[4] ] == {{(-I)*Sqrt[2]*GF*(2*gvdd - gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[3] , -S[4] , V[4] ] == {{I*Sqrt[2]*GF*(gvdd - 2*gvss + gvuu)*Vud, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[2] , V[4] ] == {{-2*fpi*GF*(gvdd - gvuu)*Vud, 0}, {0, 0}},

C[ -F[1] , F[2] , S[1] , -S[2] , V[4] ] == {{(2*I)*GF*(gvdd - gvuu)*Vud, 0}, {0, 0}},

C[ -F[1] , F[2] , S[5] , -S[4] , V[4] ] == {{I*Sqrt[3]*GF*(gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[4] , V[4] ] == {{-2*fpi*GF*(gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[4] , S[1] , V[4] ] == {{I*GF*(gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , S[3] , -S[2] , V[4] ] == {{(-I)*Sqrt[2]*GF*(2*gvdd - gvss - gvuu)*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , S[3] , S[4] , V[1] ] == {{I*Sqrt[2]*GF*qe*Vud, 0}, {0, 0}},

C[ -F[2] , F[1] , S[2] , V[1] ] == {{-2*fpi*GF*qe*Vud, 0}, {0, 0}},

C[ -F[2] , F[1] , S[1] , S[2] , V[1] ] == {{(-2*I)*GF*qe*Vud, 0}, {0, 0}},

C[ -F[2] , F[1] , S[5] , S[4] , V[1] ] == {{(-I)*Sqrt[3]*GF*qe*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , V[1] ] == {{-2*fpi*GF*qe*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , S[1] , V[1] ] == {{(-I)*GF*qe*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , -S[3] , S[2] , V[1] ] == {{I*Sqrt[2]*GF*qe*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , S[5] , S[3] , S[4] ] == {{(-2*Sqrt[2/3]*GF*Vud)/fpi, 0}, {(Sqrt[2/3]*GF*Vud)/fpi, 0}, {(Sqrt[2/3]*GF*Vud)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[5] , S[4] ] == {{(-I)*Sqrt[3]*GF*Vus, 0}, {I*Sqrt[3]*GF*Vus, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[5] , S[5] , S[4] ] == {{(GF*Vus)/fpi, 0}, {(GF*Vus)/fpi, 0}, {(-2*GF*Vus)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[5] , S[4] , S[1] ] == {{(GF*Vus)/(Sqrt[3]*fpi), 0}, {(-2*GF*Vus)/(Sqrt[3]*fpi), 0}, {(GF*Vus)/(Sqrt[3]*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[5] , -S[3] , S[2] ] == {{(Sqrt[2/3]*GF*Vus)/fpi, 0}, {(-2*Sqrt[2/3]*GF*Vus)/fpi, 0}, {(Sqrt[2/3]*GF*Vus)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[3] , S[4] ] == {{I*Sqrt[2]*GF*Vud, 0}, {(-I)*Sqrt[2]*GF*Vud, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[3] , S[4] , S[1] ] == {{-((Sqrt[2]*GF*Vud)/fpi), 0}, {(Sqrt[2]*GF*Vud)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[3] , -S[3] , S[2] ] == {{(-2*GF*Vud)/(3*fpi), 0}, {(4*GF*Vud)/(3*fpi), 0}, {(-2*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[3] , -S[3] , S[4] ] == {{(4*GF*Vus)/(3*fpi), 0}, {(-2*GF*Vus)/(3*fpi), 0}, {(-2*GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , -S[3] , S[2] ] == {{I*Sqrt[2]*GF*Vus, 0}, {(-I)*Sqrt[2]*GF*Vus, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , -S[4] , S[2] ] == {{(-2*GF*Vud)/(3*fpi), 0}, {(4*GF*Vud)/(3*fpi), 0}, {(-2*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] ] == {{2*fpi*GF*Vus, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , S[4] , -S[4] ] == {{(-4*GF*Vus)/(3*fpi), 0}, {(-4*GF*Vus)/(3*fpi), 0}, {(8*GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , S[1] ] == {{I*GF*Vus, 0}, {(-I)*GF*Vus, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , S[1] , S[1] ] == {{(-2*GF*Vus)/(3*fpi), 0}, {(GF*Vus)/(3*fpi), 0}, {(GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[4] , S[2] , -S[2] ] == {{(-2*GF*Vus)/(3*fpi), 0}, {(-2*GF*Vus)/(3*fpi), 0}, {(4*GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[1] , S[2] ] == {{(-2*I)*GF*Vud, 0}, {(2*I)*GF*Vud, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[1] , S[1] , S[2] ] == {{(4*GF*Vud)/(3*fpi), 0}, {(4*GF*Vud)/(3*fpi), 0}, {(-8*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , -S[3] , S[1] , S[2] ] == {{0, 0}, {-((Sqrt[2]*GF*Vus)/fpi), 0}, {(Sqrt[2]*GF*Vus)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[2] , F[1] , S[2] ] == {{2*fpi*GF*Vud, 0}, {0, 0}},

C[ -F[2] , F[1] , S[2] , S[2] , -S[2] ] == {{(-4*GF*Vud)/(3*fpi), 0}, {(-4*GF*Vud)/(3*fpi), 0}, {(8*GF*Vud)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[3] , -S[4] , V[1] ] == {{I*Sqrt[2]*GF*qe*Vud, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[2] , V[1] ] == {{2*fpi*GF*qe*Vud, 0}, {0, 0}},

C[ -F[1] , F[2] , S[1] , -S[2] , V[1] ] == {{(-2*I)*GF*qe*Vud, 0}, {0, 0}},

C[ -F[1] , F[2] , S[5] , -S[4] , V[1] ] == {{(-I)*Sqrt[3]*GF*qe*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[4] , V[1] ] == {{2*fpi*GF*qe*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[4] , S[1] , V[1] ] == {{(-I)*GF*qe*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , S[3] , -S[2] , V[1] ] == {{I*Sqrt[2]*GF*qe*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , S[5] , -S[3] , -S[4] ] == {{(-2*Sqrt[2/3]*GF*Vud)/fpi, 0}, {(Sqrt[2/3]*GF*Vud)/fpi, 0}, {(Sqrt[2/3]*GF*Vud)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[5] , -S[4] ] == {{I*Sqrt[3]*GF*Vus, 0}, {(-I)*Sqrt[3]*GF*Vus, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[5] , S[5] , -S[4] ] == {{(GF*Vus)/fpi, 0}, {(GF*Vus)/fpi, 0}, {(-2*GF*Vus)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[5] , -S[4] , S[1] ] == {{(GF*Vus)/(Sqrt[3]*fpi), 0}, {(-2*GF*Vus)/(Sqrt[3]*fpi), 0}, {(GF*Vus)/(Sqrt[3]*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[3] , -S[3] , -S[4] ] == {{(-2*GF*Vus)/(3*fpi), 0}, {(4*GF*Vus)/(3*fpi), 0}, {(-2*GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[3] , -S[4] ] == {{(-I)*Sqrt[2]*GF*Vud, 0}, {I*Sqrt[2]*GF*Vud, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[3] , -S[4] , S[1] ] == {{-((Sqrt[2]*GF*Vud)/fpi), 0}, {(Sqrt[2]*GF*Vud)/fpi, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , S[4] , -S[4] , -S[4] ] == {{(8*GF*Vus)/(3*fpi), 0}, {(-4*GF*Vus)/(3*fpi), 0}, {(-4*GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[4] ] == {{2*fpi*GF*Vus, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[4] , S[1] ] == {{(-I)*GF*Vus, 0}, {I*GF*Vus, 0}, {0, 0}, {0, 0}},

C[ -F[1] , F[2] , -S[4] , S[1] , S[1] ] == {{(-2*GF*Vus)/(3*fpi), 0}, {(GF*Vus)/(3*fpi), 0}, {(GF*Vus)/(3*fpi), 0}, {0, 0}, {0, 0}, {0, 0}}

}

