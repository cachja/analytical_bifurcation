(* ::Package:: *)

(* ::Package:: *)
(**)


(* Cylinder, bounded annulus: frozen Stokes base flow + 2 angular modes.
   This script computes M, D, C and the reduced matrix B[Re] from scratch. *)

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

LapPolar[f_] := D[f, {r, 2}] + D[f, r]/r + D[f, {th, 2}]/r^2;
JacPolar[a_, b_] := (D[a, r] D[b, th] - D[a, th] D[b, r])/r;

DeltaR = 1 - R^2 + (R^2 + 1) Log[R];
A0 = (1 - R^2)/(2 DeltaR);
B0 = R^2/(2 DeltaR);
C0 = (R^2 + 1)/DeltaR;
D0 = -1/(2 DeltaR);
fr[r_] := A0 r + B0/r + C0 r Log[r] + D0 r^3;
psi0[r_, th_] := fr[r] Sin[th];

ur0[r_, th_] := (1/r) D[psi0[r, th], th] // Simplify;
uth0[r_, th_] := -D[psi0[r, th], r] // Simplify;

Print["Stokes base streamfunction psi0(r,th) = f(r) sin(th)"];
Print[TraditionalForm[psi0[r, th]]];
Print["Velocity components:"];
Print[TraditionalForm[ur0[r, th]]];
Print[TraditionalForm[uth0[r, th]]];

q[r_] := (r - 1)^2 (R - r)^2;
phi[1][r_, th_] := q[r] Cos[th];
phi[2][r_, th_] := q[r] Cos[2 th];
lapPhi[j_][r_, th_] := Evaluate[LapPolar[phi[j][r, th]] // FullSimplify];

Mentry[i_, j_] := safeNIntegrate[
  ((D[phi[i][r, th], r] D[phi[j][r, th], r]) + (D[phi[i][r, th], th] D[phi[j][r, th], th])/r^2) r,
  {r, 1, R}, {th, 0, 2 Pi}
];
Dentry[i_, j_] := safeNIntegrate[
  lapPhi[i][r, th] lapPhi[j][r, th] r,
  {r, 1, R}, {th, 0, 2 Pi}
];
Centry[i_, j_] := safeNIntegrate[
  phi[i][r, th] (JacPolar[psi0[r, th], lapPhi[j][r, th]] + JacPolar[phi[j][r, th], LapPolar[psi0[r, th]]]) r,
  {r, 1, R}, {th, 0, 2 Pi}
];

M = Table[Mentry[i, j], {i, 1, 2}, {j, 1, 2}] // Chop;
Dmat = Table[Dentry[i, j], {i, 1, 2}, {j, 1, 2}] // Chop;
Cmat = Table[Centry[i, j], {i, 1, 2}, {j, 1, 2}] // Chop;
Bvisc = LinearSolve[M, Dmat] // Chop;
Bconv = -LinearSolve[M, Cmat] // Chop;
B[Re_] := (Bconv - Bvisc/Re) // Chop;

Print["M = ", MatrixForm[M]];
Print["D = ", MatrixForm[Dmat]];
Print["C = ", MatrixForm[Cmat]];
Print["Bconv = ", MatrixForm[Bconv]];
Print["Bvisc = ", MatrixForm[Bvisc]];

TrackedPair[Re_?NumericQ] := Module[{ev = Eigenvalues[N[B[Re], 30]], pos},
  pos = Select[ev, Im[#] > 10^-10 &];
  If[pos =!= {}, First@SortBy[pos, -Im[#] &],
    Module[{re = Sort[Re /@ ev], k},
      k = Ordering[Abs[Differences[re]], 1][[1]];
      (re[[k]] + re[[k + 1]])/2
    ]
  ]
];

reGrid = N@Subdivide[0.05, 120, 1200];
branchData = Table[{Re, TrackedPair[Re]}, {Re, reGrid}];
firstComplex = SelectFirst[branchData, Abs[Im[#[[2]]]] > 10^-8 &][[1]];
Print["First complex pair appears near Re = ", N[firstComplex, 8]];

Export["cylinder_stokes_2modes_branch.csv", ({#[[1]], Re[#[[2]]], Im[#[[2]]]} & /@ branchData), "CSV"];
Print["Exported branch data to cylinder_stokes_2modes_branch.csv"];

(* ===== Exact 2x2 eigenvalue formulas ===== *)
tr[Re_] := Tr[B[Re]] // FullSimplify;
detB[Re_] := Det[B[Re]] // FullSimplify;
disc[Re_] := tr[Re]^2 - 4 detB[Re] // FullSimplify;

lambdaPlus[Re_] := (tr[Re] + Sqrt[disc[Re]])/2 // FullSimplify;
lambdaMinus[Re_] := (tr[Re] - Sqrt[disc[Re]])/2 // FullSimplify;

Plot[
  Evaluate[
    {
      Re[lambdaPlus[r]],
      Re[lambdaMinus[r]],
      Im[lambdaPlus[r]],
      -Im[lambdaPlus[r]]
    }
  ],
  {r, 0, 120},
  PlotStyle -> {
    {Blue, Thick},
    {Blue, Thick},
    {Red, Dashed, Thick},
    {Red, Dashed, Thick}
  },
  AxesLabel -> {"Re", "eigenvalue components"},
  GridLines -> {{0}, {0}},
  ImageSize -> Large,
  PlotRange->{{0,120},{-0.5,0.5}}
]




