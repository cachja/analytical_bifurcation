(* ::Package:: *)

(* ::Package:: *)
(* Cylinder, bounded annulus: first-order corrected Stokes base flow + 3 angular modes. *)


ClearAll["Global`*"];

safeNIntegrate[expr_, vars__] :=
  Quiet[
    Chop[
      NIntegrate[
        expr, vars,
        Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 10000},
        AccuracyGoal -> 6,
        PrecisionGoal -> 6,
        MaxRecursion -> 8
      ],
      10^-10
    ],
    {NIntegrate::slwcon, NIntegrate::eincr}
  ];


$Assumptions = 1 < r < R && 0 <= th <= 2 Pi;

R = 5.0;
acc = 10;
prec = 10;

LapPolar[f_] := D[f, {r, 2}] + D[f, r]/r + D[f, {th, 2}]/r^2;
JacPolar[a_, b_] := (D[a, r] D[b, th] - D[a, th] D[b, r])/r;

DeltaR = 1 - R^2 + (R^2 + 1) Log[R];
A0 = (1 - R^2)/(2 DeltaR);
B0 = R^2/(2 DeltaR);
C0c = (R^2 + 1)/DeltaR;
D0 = -1/(2 DeltaR);
fr[r_] := A0 r + B0/r + C0c r Log[r] + D0 r^3;
psi0[r_, th_] := fr[r] Sin[th];
L1f[r_] := D[fr[r], {r, 2}] + D[fr[r], r]/r - fr[r]/r^2 // FullSimplify;

forcing[r_] := ((D[fr[r], r] L1f[r] - fr[r] D[L1f[r], r])/(2 r)) // FullSimplify;

geq = Derivative[4][g][r] + (2/r) Derivative[3][g][r] - (9/r^2) Derivative[2][g][r] + (9/r^3) Derivative[1][g][r] == forcing[r];
gbc = {g[1] == 0, Derivative[1][g][1] == 0, g[R] == 0, Derivative[1][g][R] == 0};

gSol = NDSolveValue[{geq, gbc}, g, {r, 1, R},
  AccuracyGoal -> acc, PrecisionGoal -> prec,
  Method -> {"BoundaryValues" -> {"Shooting"}}
];
psi1[r_, th_] := gSol[r] Sin[2 th];

Print["Solved first correction g(r). Sample values:"];
Print[Table[{rr, gSol[rr]}, {rr, 1, R, (R - 1)/4.0}]];

q[r_] := (r - 1)^2 (R - r)^2;
phi[1][r_, th_] := q[r] Cos[th];
phi[2][r_, th_] := q[r] Cos[2 th];
phi[3][r_, th_] := q[r] Cos[3 th];
lapPhi[j_][r_, th_] := Evaluate[LapPolar[phi[j][r, th]] // FullSimplify];

Mentry[i_, j_] := safeNIntegrate[
  ((D[phi[i][r, th], r] D[phi[j][r, th], r]) + (D[phi[i][r, th], th] D[phi[j][r, th], th])/r^2) r,
  {r, 1, R}, {th, 0, 2 Pi}
];
Dentry[i_, j_] := safeNIntegrate[
  lapPhi[i][r, th] lapPhi[j][r, th] r,
  {r, 1, R}, {th, 0, 2 Pi}
];
C0entry[i_, j_] := safeNIntegrate[
  phi[i][r, th] (JacPolar[psi0[r, th], lapPhi[j][r, th]] + JacPolar[phi[j][r, th], LapPolar[psi0[r, th]]]) r,
  {r, 1, R}, {th, 0, 2 Pi}
];
C1entry[i_, j_] := safeNIntegrate[
  phi[i][r, th] (JacPolar[psi1[r, th], lapPhi[j][r, th]] + JacPolar[phi[j][r, th], LapPolar[psi1[r, th]]]) r,
  {r, 1, R}, {th, 0, 2 Pi}
];

M = Table[Mentry[i, j], {i, 1, 3}, {j, 1, 3}] // Chop;
Dmat = Table[Dentry[i, j], {i, 1, 3}, {j, 1, 3}] // Chop;
C0mat = Table[C0entry[i, j], {i, 1, 3}, {j, 1, 3}] // Chop;
C1mat = Table[C1entry[i, j], {i, 1, 3}, {j, 1, 3}] // Chop;
Bvisc = LinearSolve[M, Dmat] // Chop;
Bconv0 = -LinearSolve[M, C0mat] // Chop;
Bconv1 = -LinearSolve[M, C1mat] // Chop;
B[Re_] := (Bconv0 + Re Bconv1 - Bvisc/Re) // Chop;

Print["Bconv0 = ", MatrixForm[Bconv0]];
Print["Bconv1 = ", MatrixForm[Bconv1]];
Print["Bvisc  = ", MatrixForm[Bvisc]];

PositiveImagEigen[Re_?NumericQ] := Module[{ev = Eigenvalues[N[B[Re], 40]], pos},
  pos = Select[ev, Im[#] > 10^-10 &];
  If[pos === {}, Missing["NoComplex"], First@SortBy[pos, -Im[#] &]]
];

reGrid1 = N@Subdivide[0.05, 120, 2400];
vals1 = Table[{Re, PositiveImagEigen[Re]}, {Re, reGrid1}];
firstComplex = SelectFirst[vals1, #[[2]] =!= Missing["NoComplex"] &][[1]];
Print["Complex branch appears near Re = ", N[firstComplex, 8]];
branchData = Table[
  Module[{lam = PositiveImagEigen[Re]},
    If[lam === Missing["NoComplex"], {Re, Indeterminate, 0}, {Re, Re[lam], Im[lam]}]
  ],
  {Re, reGrid1}
];
withIm = Select[branchData, NumericQ[#[[2]]] &];
crossPair = Select[Partition[withIm, 2, 1], #[[1, 2]] <= 0 && #[[2, 2]] >= 0 &];
If[crossPair =!= {}, Print["First crossing bracket ~ ", crossPair[[1, 1, 1]], "  --  ", crossPair[[1, 2, 1]]]];
Export["cylinder_corrected_3modes_branch.csv", branchData, "CSV"];
Print["Exported branch data to cylinder_corrected_3modes_branch.csv"];




(* ===== Exact 3x3 corrected cylinder matrix, find roots of cubic polynomials and plot complex eigpair ===== *)

Clear[lam, BExact, polyExact, roots3, lambdaStar, rePart, imPart];

BExact[rey_] := Rationalize[Bconv0, 0] + rey Rationalize[Bconv1, 0] - Rationalize[Bvisc, 0]/rey;

polyExact[rey_, lam_] := Expand[CharacteristicPolynomial[BExact[rey], lam]];

roots3[rey_?NumericQ] := Module[{rr = Rationalize[rey, 0]},
  N[
    Table[
      Root[Function[lam, polyExact[rr, lam]], k],
      {k, 1, 3}
    ],
    40
  ]
];

lambdaStar[rey_?NumericQ] := Module[{ev, pos},
  ev = roots3[rey];
  pos = Select[ev, Im[#] > 10^-12 &];
  If[pos === {},
    Missing["NoComplex"],
    First[pos]
  ]
];

rePart[rey_?NumericQ] := Module[{z = lambdaStar[rey]},
  If[z === Missing["NoComplex"], Indeterminate, Re[z]]
];

imPart[rey_?NumericQ] := Module[{z = lambdaStar[rey]},
  If[z === Missing["NoComplex"], 0, Im[z]]
];
Plot[
  Evaluate[{rePart[r], imPart[r], -imPart[r]}],
  {r, 0, 120},
  PlotStyle -> {
    {Blue, Thick},
    {Red, Dashed, Thick},
    {Red, Dashed, Thick}
  },
  AxesLabel -> {"Re", "eigenvalue components"},
  GridLines -> {{0.2012, 70.57}, {0}},
  PlotRange -> {{0,120},{-1,1}},
  ImageSize -> Large
]

