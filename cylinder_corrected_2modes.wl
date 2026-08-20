ClearAll["Global`*"];

(* ============================================================ *)
(* Cylinder model with a selectable second radial perturbation profile    *)
(* ============================================================ *)

(* "standard": q(r) = (r - 1)(r-R)
   "enriched": q(r) = r^2(r - 1)(r-R)
*)

basisChoice = "standard";

(* ============================================================ *)
(* Parameters *)
(* ============================================================ *)

R0 = 5.0;

reyMin = 0.05;
reyMax = 120.0;
nGrid = 2400;

accGoal = 10;
precGoal = 10;
maxRec = 14;

$Assumptions = 1 < r < R0 && Element[r, Reals];


Print["Two-mode cylinder model: q cos(theta), r^2 q cos(2 theta)"];
Print["R0 = ", R0];

(* ============================================================ *)
(* Radial operators and angular selection rules *)
(* ============================================================ *)

ClearAll[Lrad, Icc, Iss];

Lrad[n_, h_] := D[h, {r, 2}] + D[h, r]/r - (n^2/r^2) h;

Icc[m_, k_, n_] :=
  (Pi/2) ((1 - KroneckerDelta[n, k]) KroneckerDelta[m, Abs[n - k]]
    + KroneckerDelta[m, n + k]);

Iss[m_, k_, n_] :=
  (Pi/2) ((1 - KroneckerDelta[n, k]) KroneckerDelta[m, Abs[n - k]]
    - KroneckerDelta[m, n + k]);

(* ============================================================ *)
(* Stokes base psi0 = f(r) sin(theta) *)
(* ============================================================ *)

DeltaR = 1 - R0^2 + (R0^2 + 1) Log[R0];

A0 = (1 - R0^2)/(2 DeltaR);
B0 = R0^2/(2 DeltaR);
C0c = (R0^2 + 1)/DeltaR;
D0 = -1/(2 DeltaR);

fr = A0 r + B0/r + C0c r Log[r] + D0 r^3;

L1f = FullSimplify[Lrad[1, fr]];
L1fp = FullSimplify[D[L1f, r]];

forcing = FullSimplify[(D[fr, r] L1f - fr D[L1f, r])/(2 r)];

Print["Check L1 f:"];
Print[FullSimplify[L1f]];

(* ============================================================ *)
(* Corrected first base deformation psi1 = g(r) sin(2 theta) *)
(* Choice A: old J definition, so L2^2 g = -forcing *)
(* ============================================================ *)

gp = FullSimplify[
  r^6/(192 DeltaR^2)
  - ((R0^2 + 1) r^4 Log[r])/(12 DeltaR^2)
  + ((R0^4 - 6 R0^2 - 3) r^2 Log[r])/(32 DeltaR^2)
  - ((R0^2 + 1)^2 r^2 Log[r]^2)/(16 DeltaR^2)
];

gh = c1 + c2 r^2 + c3 r^4 + c4/r^2;

gGeneral = gh + gp;

bcVec = FullSimplify[
  {
    gGeneral /. r -> 1,
    D[gGeneral, r] /. r -> 1,
    gGeneral /. r -> R0,
    D[gGeneral, r] /. r -> R0
  }
];

ca = CoefficientArrays[bcVec, {c1, c2, c3, c4}];

constPart = Normal[ca[[1]]];
coeffMat = Normal[ca[[2]]];

solVec = LinearSolve[coeffMat, -constPart];

solConst = Thread[{c1, c2, c3, c4} -> solVec];

g = FullSimplify[gGeneral /. solConst];

(* Corrected sign *)
gCorr = FullSimplify[-g];

L2g = FullSimplify[Lrad[2, gCorr]];
L2gp = FullSimplify[D[L2g, r]];

Print["Check corrected g boundary conditions:"];
Print[
  Chop[
    N[
      {
        gCorr /. r -> 1,
        D[gCorr, r] /. r -> 1,
        gCorr /. r -> R0,
        D[gCorr, r] /. r -> R0
      },
      16
    ]
  ]
];

Print["Corrected g sample values:"];
Print[
  N[
    Table[{rr, gCorr /. r -> rr}, {rr, 1, R0, (R0 - 1)/4.0}],
    12
  ]
];

(* ============================================================ *)
(* Basis choice *)
(* ============================================================ *)

(* Available choices:
   "standard" : q cos(theta), q cos(2 theta)
   "enriched" : q cos(theta), r^2 q cos(2 theta)
*)

qBase = (r - 1)^2 (R0 - r)^2;

basis = Switch[
  basisChoice,
  
  "standard",
  {
    {1, 0},
    {2, 0}
  },
  
  "enriched",
  {
    {1, 0},
    {2, 2}
  },
  
  _,
  Print[
    "Unknown basisChoice: ", basisChoice,
    ". Use \"standard\" or \"enriched\"."
  ];
  Abort[]
];

nb = Length[basis];

nOf[i_] := basis[[i, 1]];
pOf[i_] := basis[[i, 2]];

h[i_] := h[i] = FullSimplify[r^pOf[i] qBase];

Lh[i_] := Lh[i] = FullSimplify[Lrad[nOf[i], h[i]]];

Lhp[i_] := Lhp[i] = FullSimplify[D[Lh[i], r]];

Print["Basis choice = ", basisChoice];
Print["Basis specification {angular level, radial power} = ", basis];

If[
  basisChoice == "standard",
  Print["phi_1 = q(r) cos(theta)"];
  Print["phi_2 = q(r) cos(2 theta)"];
];

If[
  basisChoice == "enriched",
  Print["phi_1 = q(r) cos(theta)"];
  Print["phi_2 = r^2 q(r) cos(2 theta)"];
];

(* ============================================================ *)
(* Radial integration *)
(* ============================================================ *)

ClearAll[RadialInt];

RadialInt[expr_] := RadialInt[expr] =
  Quiet[
    Chop[
      NIntegrate[
        Evaluate[expr],
        {r, 1, R0},
        AccuracyGoal -> accGoal,
        PrecisionGoal -> precGoal,
        MaxRecursion -> maxRec,
        Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
      ],
      10^-11
    ],
    {
      NIntegrate::slwcon,
      NIntegrate::eincr,
      NIntegrate::ncvb
    }
  ];

(* ============================================================ *)
(* Matrix entries *)
(* ============================================================ *)

ClearAll[MEntry, DEntry, CBaseEntry];

MEntry[i_, j_] := MEntry[i, j] =
  If[
    nOf[i] =!= nOf[j],
    0,
    Pi RadialInt[
      r D[h[i], r] D[h[j], r]
      + (nOf[i] nOf[j]/r) h[i] h[j]
    ]
  ];

DEntry[i_, j_] := DEntry[i, j] =
  If[
    nOf[i] =!= nOf[j],
    0,
    Pi RadialInt[
      r Lh[i] Lh[j]
    ]
  ];

CBaseEntry[i_, j_, k_, a_, A_, Ap_] := CBaseEntry[i, j, k] =
  Module[
    {m, n, icc, iss, radial},
    
    m = nOf[i];
    n = nOf[j];
    
    icc = Icc[m, k, n];
    iss = Iss[m, k, n];
    
    If[
      icc == 0 && iss == 0,
      0,
      radial =
        h[i] (
          (-n D[a, r] Lh[j] + n h[j] Ap) iss
          + (-k a Lhp[j] + k D[h[j], r] A) icc
        );
      RadialInt[radial]
    ]
  ];

(* ============================================================ *)
(* Assemble matrices *)
(* ============================================================ *)

Print["Assembling M..."];
M = Chop[N[Table[MEntry[i, j], {i, nb}, {j, nb}]]];

Print["Assembling D..."];
Dmat = Chop[N[Table[DEntry[i, j], {i, nb}, {j, nb}]]];

Print["Assembling C0..."];
C0mat = Chop[
  N[
    Table[
      CBaseEntry[i, j, 1, fr, L1f, L1fp],
      {i, nb},
      {j, nb}
    ]
  ]
];

Print["Assembling C1..."];
C1mat = Chop[
  N[
    Table[
      CBaseEntry[i, j, 2, gCorr, L2g, L2gp],
      {i, nb},
      {j, nb}
    ]
  ]
];

Bvisc = Chop[LinearSolve[M, Dmat]];
Bconv0 = Chop[-LinearSolve[M, C0mat]];
Bconv1 = Chop[-LinearSolve[M, C1mat]];

Print["M = "];
Print[MatrixForm[N[M, 12]]];

Print["Dmat = "];
Print[MatrixForm[N[Dmat, 12]]];

Print["C0mat = "];
Print[MatrixForm[N[C0mat, 12]]];

Print["C1mat = "];
Print[MatrixForm[N[C1mat, 12]]];

Print["Bconv0 = "];
Print[MatrixForm[N[Bconv0, 12]]];

Print["Bconv1 = "];
Print[MatrixForm[N[Bconv1, 12]]];

Print["Bvisc = "];
Print[MatrixForm[N[Bvisc, 12]]];

(* ============================================================ *)
(* Spectrum and branch data *)
(* ============================================================ *)

ClearAll[Bmat, EigVals, PositiveImagEigen];

Bmat[rey_?NumericQ] := N[Bconv0 + rey Bconv1 - Bvisc/rey];

EigVals[rey_?NumericQ] := Eigenvalues[Bmat[rey]];

PositiveImagEigen[rey_?NumericQ] := Module[
  {ev, pos},
  ev = EigVals[rey];
  pos = Select[ev, Im[#] > 10^-10 &];
  If[
    pos === {},
    Missing["NoComplex"],
    First@MaximalBy[pos, Re[#] &]
  ]
];

reGrid1 = N@Subdivide[reyMin, reyMax, nGrid];

vals1 = Table[
  {rey, PositiveImagEigen[rey]},
  {rey, reGrid1}
];

firstComplexEntry = SelectFirst[
  vals1,
  ! MatchQ[#[[2]], _Missing] &,
  Missing["NotFound"]
];

If[
  MatchQ[firstComplexEntry, _Missing],
  Print["No complex branch detected on the grid."],
  Print["Complex branch appears near Re = ", N[firstComplexEntry[[1]], 8]]
];

branchData = Table[
  Module[
    {lam = PositiveImagEigen[rey]},
    If[
      MatchQ[lam, _Missing],
      {rey, Indeterminate, 0},
      {rey, Re[lam], Im[lam]}
    ]
  ],
  {rey, reGrid1}
];

withIm = Select[branchData, NumericQ[#[[2]]] &];

crossPair = Select[
  Partition[withIm, 2, 1],
  #[[1, 2]] <= 0 && #[[2, 2]] >= 0 &
];

If[
  crossPair =!= {},
  Print[
    "First crossing bracket ~ ",
    crossPair[[1, 1, 1]],
    " -- ",
    crossPair[[1, 2, 1]]
  ]
];

outputName = Switch[
  basisChoice,
  "standard", "cylinder_2modes_q_q_branch.csv",
  "enriched", "cylinder_2modes_q_r2q_branch.csv"
];

Export[
  outputName,
  branchData,
  "CSV"
];

Print["Exported branch data to ", outputName];

(* ============================================================ *)
(* Plot tracked complex pair *)
(* ============================================================ *)
verticalGridLines = Switch[
  basisChoice,
  "standard", {20.54935},
  "enriched", {23.83}
];


ListLinePlot[
  {
    Select[branchData[[All, {1, 2}]], NumericQ[#[[2]]] &],
    Select[branchData[[All, {1, 3}]], NumericQ[#[[2]]] &],
    Select[branchData[[All, {1, 3}]], NumericQ[#[[2]]] &] /. {x_, y_} :> {x, -y}
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
  Frame -> True,
  FrameLabel -> {"Re", "eigenvalue components"},
  GridLines -> {verticalGridLines, {0}},
  PlotRange -> {{0, 35}, {-1, 1}}]
