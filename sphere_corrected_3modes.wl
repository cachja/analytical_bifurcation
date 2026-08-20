ClearAll["Global`*"];

(* ============================================================ *)
(* Sphere model with a selectable third perturbation profile    *)
(* ============================================================ *)

(* "standard": psi3(r) = (r - 1)^2/r^6
   "enriched": psi3(r) = (r - 1)^2/r^7
*)

perturbationBasisChoice = "enriched";

thirdModeExponent = Switch[
  perturbationBasisChoice,
  "standard", 6,
  "enriched", 7,
  _, Print["Unknown perturbationBasisChoice: ", perturbationBasisChoice, ". Use \"standard\" or \"enriched\"."]; Abort[]
];

Print["Perturbation basis choice = ", perturbationBasisChoice];
Print["Third-mode radial exponent = ", thirdModeExponent];

reyMin = 0.1;
reyMax = 500.0;
nGrid = 1500;

(* ============================================================ *)
(* Safe numerical integration                                   *)
(* ============================================================ *)

SetAttributes[SafeNIntegrate, HoldAll];

SafeNIntegrate[expr_, vars__] := Quiet[
  Chop[
    NIntegrate[
      Evaluate[expr], vars,
      Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 10000},
      AccuracyGoal -> 6,
      PrecisionGoal -> 6,
      MaxRecursion -> 8
    ],
    10^-10
  ],
  {NIntegrate::slwcon, NIntegrate::eincr, NIntegrate::inumr}
];

$Assumptions = r > 1 && 0 < th < Pi && 0 <= ph <= 2 Pi;

(* ============================================================ *)
(* Basic spherical vector operations                            *)
(* ============================================================ *)

DotSph[A_, B_] := Simplify[A[[1]] B[[1]] + A[[2]] B[[2]] + A[[3]] B[[3]]];

CurlSph[A_] := Simplify[
  {
    (D[Sin[th] A[[3]], th] - D[A[[2]], ph])/(r Sin[th]),
    (D[A[[1]], ph]/Sin[th] - D[r A[[3]], r])/r,
    (D[r A[[2]], r] - D[A[[1]], th])/r
  }
];

DirDerSph[V_, A_] := Module[{Vr, Vt, Vp, Ar, At, Ap, T},
  Vr = V[[1]];
  Vt = V[[2]];
  Vp = V[[3]];
  Ar = A[[1]];
  At = A[[2]];
  Ap = A[[3]];
  T[f_] := Vr D[f, r] + (Vt/r) D[f, th] + (Vp/(r Sin[th])) D[f, ph];
  Simplify[
    {
      T[Ar] - (Vt At + Vp Ap)/r,
      T[At] + (Vt Ar - Vp Ap Cot[th])/r,
      T[Ap] + (Vp Ar + Vp At Cot[th])/r
    }
  ]
];

(* ============================================================ *)
(* Frozen Stokes base flow                                      *)
(* ============================================================ *)

urS[r_, th_] := Cos[th] (1 - 3/(2 r) + 1/(2 r^3));
uthS[r_, th_] := -Sin[th] (1 - 3/(4 r) - 1/(4 r^3));
uphS[r_, th_] := 0;

US = {urS[r, th], uthS[r, th], uphS[r, th]};

Print["Frozen Stokes base flow U_S in spherical components:"];
Print[TraditionalForm[US]];

Fforce = FullSimplify[DirDerSph[US, US]];

Print["Forcing F = (U_S . grad) U_S in spherical components:"];
Print[TraditionalForm[Fforce]];

(* ============================================================ *)
(* General poloidal mode formula                                *)
(* ============================================================ *)

PoloidalMode[l_Integer, psiFun_, Y_] := Module[{alpha},
  alpha = Simplify[D[r psiFun[r], r]/r];
  Simplify[
    {
      l (l + 1) psiFun[r] Y/r,
      alpha D[Y, th],
      alpha D[Y, ph]/Sin[th]
    }
  ]
];

(* ============================================================ *)
(* PART A: reduced first-order base correction                  *)
(* ============================================================ *)

Yb1[th_] := Cos[th];
Yb2[th_] := (3 Cos[th]^2 - 1)/2;

psib1[r_] := (r - 1)^2/r^4;
psib2[r_] := (r - 1)^2/r^5;

BaseMode1 = FullSimplify[PoloidalMode[1, psib1, Yb1[th]]];
BaseMode2 = FullSimplify[PoloidalMode[2, psib2, Yb2[th]]];

CurlBaseMode1 = FullSimplify[CurlSph[BaseMode1]];
CurlBaseMode2 = FullSimplify[CurlSph[BaseMode2]];

baseBasis = {BaseMode1, BaseMode2};
curlBaseBasis = {CurlBaseMode1, CurlBaseMode2};

Print["Axisymmetric correction basis BaseMode1 = "];
Print[TraditionalForm[BaseMode1]];

Print["Axisymmetric correction basis BaseMode2 = "];
Print[TraditionalForm[BaseMode2]];

vol = r^2 Sin[th];

PhiIntegrated[expr_] := FullSimplify[
  Integrate[expr, {ph, 0, 2 Pi}],
  Assumptions -> r > 1 && 0 < th < Pi
];

Sentry[i_, j_] := Module[{integrand},
  integrand = PhiIntegrated[DotSph[curlBaseBasis[[j]], curlBaseBasis[[i]]] vol];
  SafeNIntegrate[integrand, {r, 1, Infinity}, {th, 0, Pi}]
];

RHSentry[i_] := Module[{integrand},
  integrand = PhiIntegrated[-DotSph[Fforce, baseBasis[[i]]] vol];
  SafeNIntegrate[integrand, {r, 1, Infinity}, {th, 0, Pi}]
];

Smat = Chop[Table[Sentry[i, j], {i, 1, 2}, {j, 1, 2}]];
RHS = Chop[Table[RHSentry[i], {i, 1, 2}]];
coeff = Chop[LinearSolve[Smat, RHS]];

Print["Reduced Stokes correction matrix S = "];
Print[MatrixForm[Smat]];

Print["Reduced Stokes forcing RHS = "];
Print[MatrixForm[RHS]];

Print["Correction coefficients {c1,c2} = ", coeff];

U1red = FullSimplify[coeff[[1]] BaseMode1 + coeff[[2]] BaseMode2];

Print["Reduced first correction U1red = "];
Print[TraditionalForm[U1red]];

(* ============================================================ *)
(* PART B: three-mode real m=1 perturbation chain               *)
(* ============================================================ *)

Y1[th_, ph_] := Sin[th] Cos[ph];
Y2[th_, ph_] := Sin[th] Cos[th] Cos[ph];
Y3[th_, ph_] := Sin[th] (5 Cos[th]^2 - 1) Cos[ph];

psi1[r_] := (r - 1)^2/r^4;
psi2[r_] := (r - 1)^2/r^5;
psi3[r_] := (r - 1)^2/r^thirdModeExponent;

Print["Perturbation radial profiles:"];
Print["psi1(r) = (r - 1)^2/r^4"];
Print["psi2(r) = (r - 1)^2/r^5"];
Print["psi3(r) = (r - 1)^2/r^", thirdModeExponent];

boundaryCheck = Chop[
  N[
    {
      psi1[1], Derivative[1][psi1][1],
      psi2[1], Derivative[1][psi2][1],
      psi3[1], Derivative[1][psi3][1]
    },
    16
  ]
];

Print["Boundary-condition check = ", boundaryCheck];

Phi1 = FullSimplify[PoloidalMode[1, psi1, Y1[th, ph]]];
Phi2 = FullSimplify[PoloidalMode[2, psi2, Y2[th, ph]]];
Phi3 = FullSimplify[PoloidalMode[3, psi3, Y3[th, ph]]];

CurlPhi1 = FullSimplify[CurlSph[Phi1]];
CurlPhi2 = FullSimplify[CurlSph[Phi2]];
CurlPhi3 = FullSimplify[CurlSph[Phi3]];

basis = {Phi1, Phi2, Phi3};
curlBasis = {CurlPhi1, CurlPhi2, CurlPhi3};
nb = Length[basis];

Print["Phi1 = "];
Print[TraditionalForm[Phi1]];

Print["Phi2 = "];
Print[TraditionalForm[Phi2]];

Print["Phi3 = "];
Print[TraditionalForm[Phi3]];

(* ============================================================ *)
(* Reduced matrix entries                                       *)
(* ============================================================ *)

Mentry[i_, j_] := Module[{integrand},
  integrand = PhiIntegrated[DotSph[basis[[j]], basis[[i]]] vol];
  SafeNIntegrate[integrand, {r, 1, Infinity}, {th, 0, Pi}]
];

Dentry[i_, j_] := Module[{integrand},
  integrand = PhiIntegrated[DotSph[curlBasis[[j]], curlBasis[[i]]] vol];
  SafeNIntegrate[integrand, {r, 1, Infinity}, {th, 0, Pi}]
];

K0entry[i_, j_] := Module[{convShear, integrand},
  convShear = FullSimplify[-(DirDerSph[US, basis[[j]]] + DirDerSph[basis[[j]], US])];
  integrand = PhiIntegrated[DotSph[convShear, basis[[i]]] vol];
  SafeNIntegrate[integrand, {r, 1, Infinity}, {th, 0, Pi}]
];

K1entry[i_, j_] := Module[{convShear, integrand},
  convShear = FullSimplify[-(DirDerSph[U1red, basis[[j]]] + DirDerSph[basis[[j]], U1red])];
  integrand = PhiIntegrated[DotSph[convShear, basis[[i]]] vol];
  SafeNIntegrate[integrand, {r, 1, Infinity}, {th, 0, Pi}]
];

Print["Assembling M..."];
M = Chop[Table[Mentry[i, j], {i, 1, nb}, {j, 1, nb}]];

Print["Assembling D..."];
Dmat = Chop[Table[Dentry[i, j], {i, 1, nb}, {j, 1, nb}]];

Print["Assembling K0..."];
K0mat = Chop[Table[K0entry[i, j], {i, 1, nb}, {j, 1, nb}]];

Print["Assembling K1..."];
K1mat = Chop[Table[K1entry[i, j], {i, 1, nb}, {j, 1, nb}]];

Print["M = "];
Print[MatrixForm[M]];

Print["D = "];
Print[MatrixForm[Dmat]];

Print["K0 = "];
Print[MatrixForm[K0mat]];

Print["K1 = "];
Print[MatrixForm[K1mat]];

(* ============================================================ *)
(* Reduced spectral matrix                                      *)
(* ============================================================ *)

Bvisc = Chop[LinearSolve[M, Dmat]];
Bconv0 = Chop[LinearSolve[M, K0mat]];
Bconv1 = Chop[LinearSolve[M, K1mat]];

Bmat[rey_?NumericQ] := Chop[Bconv0 + rey Bconv1 - Bvisc/rey];

Print["Bconv0 = "];
Print[MatrixForm[Bconv0]];

Print["Bconv1 = "];
Print[MatrixForm[Bconv1]];

Print["Bvisc = "];
Print[MatrixForm[Bvisc]];

(* ============================================================ *)
(* Branch tracking                                               *)
(* ============================================================ *)

TrackedBranch[rey_?NumericQ] := Module[{ev, pos, revals, idx},
  ev = Eigenvalues[N[Bmat[rey], 40]];
  pos = Select[ev, Im[#] > 10^-10 &];

  If[
    pos =!= {},
    First@MaximalBy[pos, Re[#] &],
    revals = Sort[Re /@ ev];
    idx = Ordering[Abs[Differences[revals]], 1][[1]];
    (revals[[idx]] + revals[[idx + 1]])/2
  ]
];

reGrid = N@Subdivide[reyMin, reyMax, nGrid];

branchData = Table[
  Module[{lam},
    lam = TrackedBranch[rey];
    {rey, Re[lam], Im[lam]}
  ],
  {rey, reGrid}
];

oscCandidates = Select[branchData, Abs[#[[3]]] > 10^-8 &];

If[
  oscCandidates =!= {},
  Reosc = oscCandidates[[1, 1]];
  Print["First complex branch appears near Re = ", N[Reosc, 8]]
];

complexData = Select[branchData, Abs[#[[3]]] > 10^-8 &];

If[
  complexData =!= {},
  nearHopf = First@MaximalBy[complexData, #[[2]] &];
  Print["Closest approach to axis near Re = ", nearHopf[[1]]];
  Print["lambda(Re) approximately ", nearHopf[[2]], " + ", nearHopf[[3]], " i"]
];

(* ============================================================ *)
(* First Hopf crossing                                           *)
(* ============================================================ *)

crossingPairs = {};

If[
  Length[complexData] >= 2,
  crossingPairs = Select[
    Partition[complexData, 2, 1],
    #[[1, 2]] #[[2, 2]] < 0 &
  ]
];

If[
  crossingPairs === {},
  Print["No finite-Reynolds-number Hopf crossing detected on the scan interval."],
  firstBracket = crossingPairs[[1]];

  hopfRe = rey /. FindRoot[
    Re[TrackedBranch[rey]] == 0,
    {rey, firstBracket[[1, 1]], firstBracket[[2, 1]]}
  ];

  hopfLambda = TrackedBranch[hopfRe];
  hopfSpectrum = Eigenvalues[N[Bmat[hopfRe], 30]];

  Print["First Hopf crossing at Re = ", N[hopfRe, 12]];
  Print["Critical eigenvalue = ", N[hopfLambda, 12]];
  Print["Critical frequency omega = ", N[Im[hopfLambda], 12]];
  Print["Diameter-based Reynolds number = ", N[2 hopfRe, 12]];
  Print["Diameter-based Strouhal number omega/Pi = ", N[Im[hopfLambda]/Pi, 12]];
  Print["Full spectrum at crossing = ", N[hopfSpectrum, 12]]
];

(* ============================================================ *)
(* Export branch data                                            *)
(* ============================================================ *)

branchOutputFile = Switch[
  perturbationBasisChoice,
  "standard", "sphere_3modes_standard_branch.csv",
  "enriched", "sphere_3modes_enriched_branch.csv"
];

numericBranchData = Select[
  branchData,
  NumericQ[#[[1]]] && NumericQ[#[[2]]] && NumericQ[#[[3]]] &
];

Export[
  branchOutputFile,
  Prepend[numericBranchData, {"Re", "real", "imag"}],
  "CSV"
];

Print["Exported branch data to ", branchOutputFile];

(* ============================================================ *)
(* Plot                                                          *)
(* ============================================================ *)

plotTitle = Switch[
  perturbationBasisChoice,
  "standard", "Sphere three-mode model: standard radial basis",
  "enriched", "Sphere three-mode model: enriched third radial mode"
];

ListLinePlot[
  {
    ({#[[1]], #[[2]]} & /@ branchData),
    ({#[[1]], #[[3]]} & /@ branchData),
    ({#[[1]], -#[[3]]} & /@ branchData)
  },
  PlotStyle -> {
    {Blue, Thick},
    {Red, Dashed, Thick},
    {Red, Dashed, Thick}
  },
  PlotLegends -> {
    "Re(lambda)",
    "Im(lambda)",
    "-Im(lambda)"
  },
  PlotLabel -> plotTitle,
  AxesLabel -> {"Re", "eigenvalue components"},
  GridLines -> {
    Join[
      If[ValueQ[Reosc], {Reosc}, {}],
      If[ValueQ[hopfRe], {hopfRe}, {}]
    ],
    {0}
  },
  PlotRange -> {{0, reyMax}, {-0.5, 0.5}},
  ImageSize -> Large
]
