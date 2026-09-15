# A fast, colonising life history predicts herbicide resistance in arable weeds

Full analysis pipeline for predicting the **herbicide-resistance status of German arable weeds from
whole-plant ecological strategy** (functional traits + life history + biogeography), with phylogenetic
correction and a trait-based risk watch-list.

This repository accompanies the manuscript:

> Binama, B., Scholten, S., *et al.* *A fast, colonising life history predicts herbicide resistance in
> arable weeds.* (2026).

**Headline result.** A species' whole-plant ecological strategy — a *fast, fecund, colonising* syndrome
(annual life history, C4 photosynthesis, tall stature, high fecundity, acquisitive high-SLA leaves) —
separates resistant from non-resistant weeds well above chance (out-of-bag AUC ≈ 0.85; permutation
*P* < 0.005). The signal is robust to missing-data handling, to the learner (random forest, elastic-net and
gradient boosting agree), and to shared ancestry (phylogenetic *D*-statistic and phylogenetic logistic
regression). The applied output is a ranked watch-list of currently non-resistant weeds whose strategy most
resembles resistant ones.

---

## What this repository contains

`Reproducibility_Full_Pipeline.R` is a **single-file concatenation of the entire analysis**, assembled in
run order. Every script the manuscript's results depend on is included in one file:

- the `.R` sources are included **verbatim**;
- the two literate sources (`Herbicide_Code.qmd`, `02_Master_Analysis.Rmd`) are included as their
  **extracted R code** (via `knitr::purl`). Their fully annotated originals — with the interpretation text —
  are the `.qmd` / `.Rmd` files in the manuscript's source repository.

The file is organised into eleven stages. Stages 1–7 write intermediate outputs that Stage 8 (the main
analysis) reads; Stage 9 builds the figures and Stage 10 knits the report.

| Stage | Section | Source script | What it does |
|---|---|---|---|
| 0 | Dependencies / session setup | `install_dependencies.R` | Installs and reports on all required packages |
| 1 | Occurrence cleaning + species assembly | `Herbicide_Code.qmd` | Cleans GBIF records, assembles the 575-species table and traits |
| 2 | Geographic-range covariates | `_fix_traits.R` | EOO, AOO and altitude covariates |
| 3 | Agraratlas overlay | `07_agraratlas_overlay.R` | Agricultural-context overlay |
| 4 | Crop-data overlay | `08_cropdata_overlay.R` | District crop-share context |
| 5 | HR-status spatial surfaces | `09_hr_surface.R` | "Realized HR" surface from occurrences → map figures |
| 6 | Sites-of-action / mechanism | `10_moa_analysis.R` | HRAC mechanism composition + selection-history analysis |
| 7 | Checks | `11_checks.R` | Recording-effort control, down-sampling learning curve |
| 8 | **Main analysis** | `02_Master_Analysis.Rmd` | Modules/PCA, balanced random forest, importance, partial dependence, assessment, three-class + watch-list, and all robustness checks (phylogeny, alternative learners, SHAP) |
| 9 | Figure builders | `_make_figs.R`, `_build_fig_ms4.R`, `_compose_supp.R` | Main and supplementary figure composites |
| 10 | Render driver | `_render.R` | Knits the master report |

---

## Requirements

- **R ≥ 4.5** (developed and tested on **R 4.5.2**).
- **Pandoc** (bundled with RStudio, or install standalone) and a **LaTeX engine** — `xelatex` for the
  walkthrough, `pdflatex` for the report/supplement — only needed for Stage 10 (rendering). If you have no
  system LaTeX: `install.packages("tinytex"); tinytex::install_tinytex()`.
- **R packages** — Stage 0 of the file installs everything. To install ahead of time, run that block on its
  own (it is a self-contained copy of `install_dependencies.R`).

  > **Note:** `V.PhyloMaker2` is **not on CRAN**; Stage 0 pulls it from GitHub
  > (`jinyizju/V.PhyloMaker2`) via `remotes::install_github()`.

---

## Input data

The pipeline expects the following inputs. Place them where the code expects them, or edit the paths at the
top of the relevant stage. See the manuscript's **Data availability** statement for full citations and
access terms.

**Required for the model**

1. **`merged_data_final.csv`** — the main input: **one row per species (575 species)** with the functional
   traits, life-history and biogeographic covariates, taxonomic `family`, and the resistance response
   `HR_Status` (`HR_local` = resistant in Germany; `HR_elsewhere` = resistant abroad only;
   `Non_HR_local` = not resistant). This single file drives both the model and the phylogeny build.

**Optional (only for specific side-analyses / maps)**

2. **GBIF occurrences** — cleaned coordinates per species. Used for the spatial "realized HR" surface and
   watch-list maps. Source: GBIF, **doi:10.15468/dl.9j88hn**.
3. **Global herbicide-resistance database** — used for the mechanism-of-action analysis (Stage 6). Source:
   the International Herbicide-Resistant Weed Database (Heap, [weedscience.org](https://www.weedscience.org)).
4. **Habitat / context layers** — EVA vegetation plots, TRY traits, the BVL herbicide register,
   EUNIS/FloraVeg habitat classification, district crop statistics, and CIAT/GRIN crop-wild-relative lists,
   as described in the manuscript.

All other tables consumed by the pipeline are **regenerated** by earlier stages.

---

## How to run

`Reproducibility_Full_Pipeline.R` is a linear concatenation intended to be executed **stage by stage in
order**, because later stages read the intermediate files that earlier stages write. The recommended path:

```r
# 0. install packages (run the Stage 0 block once)
#    or, equivalently:  Rscript install_dependencies.R

# 1–7. run each stage in order; they write merged_data_final.csv, the occurrence
#      CSVs, moa_analysis.rds, reviewer_checks.rds and Figures/* used downstream

# 8.  the main analysis: models, importance, phylogeny, watch-list, nested model, CWR test
#     -> writes results_summary.rds (the single object of headline numbers)

# 9.  figure builders -> Figures/ (600-dpi PNG + vector PDF)

# 10. render driver -> Herbicide_HR_Report
```

**Graceful degradation.** The main analysis is defensive: chunks that need a cached artifact are guarded by
`eval = file.exists(...)`, and the phylogeny falls back to a taxonomy tree if the megatree is absent. Stage 8
therefore runs end-to-end even if you skip the optional map/mechanism stages — it simply omits those panels.

---

## Outputs

- `Figures/` — 600-dpi PNG + vector (`cairo_pdf`) PDF.
- `Tables/` — CSV analysis tables.
- `Models/` — fitted `.rds` objects.
- `results_summary.rds` — the single object of headline numbers consumed by the report, manuscript, deck and
  supplement.
- Exported analysis tables (including `phylo_species_data.csv`, `phylo_family_load.csv` and the Newick trees
  behind the phylogeny figures).

---

## Reproducibility notes

- **Randomness is controlled.** A global `set.seed(42)` plus local seeds cover the balanced-bootstrap forest,
  label permutations, repeated stratified cross-validation, and `missForest` imputation, so results are
  deterministic given the same package versions.
- **Run modes.** A `full` run (`ntree = 1000`, 999 permutations, leave-one-out, SHAP) reproduces the
  submission numbers; a `quick` preset (`ntree = 300`, 49 permutations, no LOO/SHAP) is a fast smoke-test
  that preserves the pipeline shape.

---

## Software manifest (tested versions)

R 4.5.2. Key package versions:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| tidyverse | 2.0.0 | | randomForest | 4.7-1.2 |
| data.table | 1.18.2.1 | | partykit | 1.2-27 |
| janitor | 2.2.1 | | glmnet | 4.1-10 |
| missForest | 1.6.1 | | xgboost | 3.2.1.1 |
| FactoMineR | 2.13 | | caret | 7.0-1 |
| factoextra | 1.0.7 | | pROC | 1.19.0.1 |
| vegan | 2.7-2 | | kernelshap | 0.9.1 |
| ape | 5.8-1 | | shapviz | 0.10.3 |
| phylolm | 2.6.5 | | kableExtra | 1.4.0 |
| caper | 1.0.4 | | patchwork | 1.3.2 |
| V.PhyloMaker2 | 0.1.0 | | sf | 1.1-1 |
| igraph | 2.2.2 | | terra | 1.9-25 |

(Full list installed by Stage 0.)

---

## Data availability & licensing

- The global herbicide-resistance records derive from the International Herbicide-Resistant Weed Database
  (Heap, weedscience.org). District crop-share statistics are official German regional statistics. Occurrence
  records are from GBIF (doi:10.15468/dl.9j88hn). Trait data are from the TRY database and vegetation data
  from EVA. Please cite the original sources when reusing the derived tables. See the manuscript for full
  methods, citations and acknowledgements.
- The **code** in this repository is released under the **MIT License** (see [`LICENSE`](LICENSE)). Input
  datasets remain under the terms of their original providers.

---

## Citation

If you use this code, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)):

> Binama, B., Scholten, S., *et al.* *A fast, colonising life history predicts herbicide resistance in
> arable weeds.* (2026).

## Contact

Blaise Binama — blaise.binama@uni-goettingen.de
Department of Crop Sciences & Centre for Integrated Breeding Research (CiBreed),
Georg-August-University Göttingen.
