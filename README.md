# TOPologically-constrained Thermodynamics-based Flux Analysis (Top-TFA)

This MATLAB code base accompanies the manuscript *"Thermodynamically feasible mass-balanced
metabolic states constrained by network topology and a linear flux-force relationship"*
(Zapararte, Marcellin, Nielsen and Saa). The repository contains the scripts,
functions, models and pre-computed data required to reproduce the two case studies of the
work and the figures derived from them.

---

## Overview

At a high level, the code implements a workflow that combines the three steps of Top-TFA using a
constraint-based metabolic model:

1. **Network-topology enumeration** — an exhaustive, loopless enumeration of all
   thermodynamically consistent reaction directionalities (flux topes) that a model
   can adopt.
2. **Thermodynamic and mass-balance constraints** — Gibbs free-energy of reaction (ΔGr),
   metabolite-concentration and charge/proton constraints layered on top of the classic
   TMFA formulation, evaluated topology by topology.
3. **A linear flux–force relationship** — an approximation of the coupling between reaction
   flux and its thermodynamic driving force, whose parameter is calibrated separately.

These ingredients are used to (i) determine which topologies are feasible, (ii) tighten the
flux and ΔGr variability ranges, (iii) sample the feasible solution space, and (iv) perform
a bottleneck / thermodynamic-sensitivity analysis. The pipeline is illustrated on two
biological systems that give the two main case-study folders.

---

## Repository layout

```
Top-TFA/
├── slopeParamDefinition.m     Calibration of the linear flux–force (slope) parameter
├── DFS/                       Topology (directionality) enumeration engine
├── looplessFxns/              Loopless flux-sampling library
├── ecoli/                     Case study 1: E. coli core model
├── clostridium/               Case study 2: Clostridium autoethanogenum model
├── manuscript.pdf, SI.pdf     Reference documents (context only)
├── LICENSE                    MIT License
└── README.md                  This file
```

### `DFS/` — Topology enumeration algorithm

Enumerates every feasible reaction directionality of a model using a
depth-first search (DFS) over the reversible reactions, pruned by loopless and
mass-balance infeasibility patterns.

- `findTopologiesNewClusters3.m` — It builds the loopless MILP structure, runs a loopless FVA,
  constructs the directionality search space, identifies infeasible balance/loop patterns,
  and compiles the full set of valid topologies (the *topology / feasibility matrix*, `tfm`).
- `fxns/` — supporting routines: loopless structure builders (`looplessStructureLP/MILP.m`),
  loopless FVA (`looplessFVA.m`), sparse null-space computation (`fastSNP.m`,
  `nullSparseBasisStructure_GUROBI.m`), the DFS pattern searches
  (`depthSearch*`), and pattern compilation / translation helpers
  (`compilePatternsRev.m`, `TranslateBinaryCfgs.m`, `readBinaryCfgs.m`, etc.).
- `C fxns/` and `juicy3_true_binary.c` — C/MEX helpers (`.c` sources with pre-compiled
  `.mexw64` binaries for Windows) that expand infeasibility patterns into the final list of
  binary configurations efficiently. Auxiliary/analysis scripts (`Testing.m`,
  `CheckingClusters.m`, `ForCytoscape.m`) help inspect and export the enumerated topologies.

### `looplessFxns/` — Flux-sampling library

A self-contained Monte-Carlo sampling library (based on loopless flux sampling,
https://github.com/SysBioengLab/looplessFluxSampler) used by both case studies to draw random
points from the feasible thermodynamic solution space. Key entry points are
`looplessFluxSampler.m` / `newLooplessFluxSampler.m`, with the samplers `ADSB.m`
(Adaptive Directions Sampling on a Box), `ll_ACHRB.m` and `HitAndRun.m`, plus problem
builders (`buildLooplessProblem.m`, `buildLinearProblem.m`), convergence diagnostics
(`psrf.m`) and utilities (`bringToBoundary.m`, `sampleSet.m`, `fast_snp.m`).

### `ecoli/` — Case study 1 (E. coli core)

Relevant code used to illustrate the capabilities on the *E. coli* core model.

- `main_e_coli_core.m` — **main entry point**. Runs the three phases end to end:
  (1) model curation (adding thermodynamic data and reducing the model into the working
  sub-models), (2) pre-processing (topology search followed by TMFA and topology-constrained
  TMFA, `Top-TMFA`, over every sub-model and dataset), and (3) sampling.
- `reducing_e_coli_core.m` — exploratory script documenting how the reduced sub-model
  (glycolysis / PPP / TCA subset with the relevant transport and exchange reactions) was
  carved out of the full core model.
- `dGr_Script.m` — computes the reaction ΔGr ranges by three increasingly constrained
  methods (concentration-only, Top-TMFA, and Top-TMFA + maximum flux) and exports the
  comparison to Excel.
- `fnxs/` — the function library for this case study: thermodynamic model construction
  (`addThermoFields.m`, `buildTherModels.m`, `newbuildTherModels.m`), model reduction
  (`modelReduction.m`, `compressModel.m`), feasibility / variability analysis
  (`checkFeasibility.m`, `FVA.m`, `MVA.m`, `newMVA.m`, `maxThermoFlux.m`), directionality
  fixing (`fixDrxns.m`, `rmvLoopyRxns.m`, `reorderTFM.m`), the sampler wrapper
  (`topTFAsampler.m`, `tcTMFAsampler.m`) and the bottleneck analysis (`bottleneckAnalysis.m`).
- `graphFxns/` — plotting scripts that reproduce the figures (`figure_5.m`, `dGrGraph.m`,
  `newDGRgraph.m`, `fluxGraph.m`, `plotDGRvsflux.m`, `segregationGraph.m`,
  `CHRRvsparallelADS.m`, `ratiovsTop.m`, and helpers `dscatter.m`, `fcsread.m`).
- `models/`, `files/`, `Results/` — input models (`e_coli_core.mat`, `COREv1/2.mat`, the
  reduced `subModels/`), intermediate `.mat` data per sub-model, and the final Excel/`.mat`
  results.

### `clostridium/` — Case study 2 (C. autoethanogenum)

Same methodology applied to a metabolic model of *Clostridium autoethanogenum*,
under experimentally measured scenarios.

- `main_clostridium.m` — **main entry point** (a function taking a measurement flag and
  sampling options). Builds the thermodynamic model, runs the topology search, discards
  non-mass-feasible topologies, runs topology-constrained TMFA (`closa_checkFeasibility.m`),
  saves the results and launches sampling.
- `mainScript.m` — convenience driver that runs repeated uniform vs. non-uniform sampling
  experiments from the saved data.
- `bottleScript.m` — runs the bottleneck / thermodynamic-sensitivity analysis for both
  measurement scenarios and writes the results to Excel.
- `fnxs/` — the `closa_*` counterparts of the E. coli functions (model construction,
  compression, feasibility checking, sampling), together with the statistical machinery for
  non-uniform / density-weighted sampling (`nonUniformADSB.m`, `metropolisHitAndRun.m`,
  `modelMVN.m`, `nearestSPD.m`, `myMVNPDF.m`, etc.).
- `graphFxns/` — figure scripts specific to this case study (`plotFig6.m`, `plotFig7.m`,
  `plotAllData.m`, `plotComparison.m`, `uniVSnonuni.m`, `corrplotScript.m`, …).
- `models/`, `files/`, `results/` — the model (`rcaModel.mat`, `thermo_rca.mat`),
  thermodynamic/concentration input data, and the pre-computed `DATA_v1/v2.mat` outputs.

### `slopeParamDefinition.m`

Stand-alone script that calibrates the parameter of the linear flux–force approximation. It
compares the exact non-linear flux–force curve against its piecewise-linear approximation
across temperatures, selects the (pseudo-)optimal slope parameter, and quantifies the
resulting approximation error.

---

## General rational

Both case studies follow the **same four-stage pipeline**, which explains the parallel
folder structure (`ecoli/` mirrors `clostridium/`, with the latter prefixing its functions
with `closa_`):

1. **Model curation** — attach thermodynamic data (ΔGf, charges, concentration bounds) and,
   for E. coli, reduce the model to a tractable working sub-network.
2. **Topology enumeration** — call the shared `DFS/` engine to enumerate all feasible
   directionalities (`tfm`).
3. **Feasibility & pre-processing** — for each topology, fix reaction directions and solve
   the TMFA / Top-TMFA problems to test feasibility and obtain tightened flux and ΔGr ranges.
4. **Sampling & analysis** — sample the feasible space with the `looplessFxns/` library and
   run the downstream ΔGr, bottleneck and figure-generation scripts.

The two leaf-level libraries (`DFS/` and `looplessFxns/`) are deliberately kept generic and
shared, while everything model-specific lives inside the corresponding case-study folder.

---

## Requirements

- **MATLAB** (R2019b or newer recommended). The following toolboxes are used across the
  scripts: Optimization Toolbox (`linprog`, `optimoptions`), Symbolic Math Toolbox
  (`slopeParamDefinition.m`), Statistics and Machine Learning Toolbox, and the Parallel
  Computing Toolbox (the samplers use `parallelFlag = 1`).
- **[COBRA Toolbox](https://opencobra.github.io/cobratoolbox/)** — for model handling and
  FVA (`initCobraToolbox`, `fluxVariability`, `findExcRxns`, …). Initialize it before
  running the main scripts (the `initCobraToolbox` calls are present but commented out).
- **[Gurobi](https://www.gurobi.com/)** — the topology-enumeration engine and null-space
  routines call Gurobi directly for the LP/MILP problems (`fastSNP.m`,
  `nullSparseBasisStructure_GUROBI.m`, the `depthSearch*` functions).
- **A C compiler / MEX** — only needed if the pre-compiled `DFS/*.mexw64` binaries (built for
  64-bit Windows) must be re-generated on another platform; recompile the `.c` sources with
  `mex`.

---

## Reproducing the key results

> Run everything from within the relevant case-study folder so that the relative
> `addpath` calls at the top of each script resolve correctly. Make sure MATLAB is started
> with COBRA and Gurobi on the path.

### Case study 1 — E. coli core

```matlab
cd ecoli
% initCobraToolbox(false)      % run once per MATLAB session
main_e_coli_core               % Parts 1–3: curation, topology + (Top-)TMFA, sampling
dGr_Script                     % ΔGr comparison across the three methods -> Results/Excels
bottleneckAnalysis             % bottleneck / thermodynamic-sensitivity analysis
```

Figures are then produced from the scripts in `ecoli/graphFxns/` (e.g. run `figure_5`,
`newDGRgraph`, `fluxGraph`). Pre-computed intermediate data is available under
`ecoli/files/` and `ecoli/Results/`, so individual downstream scripts can be run without
re-executing the full pipeline.

### Case study 2 — C. autoethanogenum

```matlab
cd clostridium
% initCobraToolbox(false)
% medition = 0/1 selects the measurement scenario; options.name labels the sampling run
options.name = "run1";
main_clostridium(1, options)   % build model, enumerate topologies, Top-TMFA, sample
mainScript                     % (optional) repeated uniform vs non-uniform sampling
bottleScript                   % bottleneck / thermodynamic-sensitivity analysis -> Excel
```

Figures are generated from `clostridium/graphFxns/` (e.g. `plotFig6`, `plotFig7`,
`uniVSnonuni`). As above, `clostridium/results/DATA_v1.mat` and `DATA_v2.mat` hold the
pre-computed outputs used by the analysis and plotting scripts.

### Flux–force parameter

```matlab
run slopeParamDefinition       % calibrates the flux–force slope parameter and its error
```

---

## Notes

- The topology enumeration is the most computationally intensive stage; on larger models it
  can take a long time and benefits from Gurobi and the compiled MEX helpers.
- Some scripts contain hard-coded absolute paths (e.g. in `reducing_e_coli_core.m`) reflecting
  the original development environment; adjust the `addpath` lines to your local layout.
- Intermediate and final results are provided in the `files/`, `models/`, `Results/` and
  `results/` subfolders so that the figures and tables can be regenerated without rerunning
  the complete pipeline.

## License

Released under the MIT License — see [LICENSE](LICENSE). © 2024 Systems Bioengineering
Laboratory.
