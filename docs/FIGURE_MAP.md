# Figure reproduction map

This map uses the figure numbering in the uploaded manuscript (`manuscript.pdf`)
and supplementary materials (`sm.pdf`). Canonical entry points are named
`manuscript_figure_<number>.R` and `supplementary_figure_<number>.R`.
Each script writes a PDF with the matching `manuscript_figure_<number>.pdf` or
`supplementary_figure_<number>.pdf` name under `figures/generated/`.

## Main manuscript

| Manuscript figure | Caption/topic | Reproduction script | Output |
| --- | --- | --- | --- |
| Figure 1 | Identifiability of optimal treatment decisions | `scripts/manuscript_figure_1.R` | `figures/generated/manuscript_figure_1.pdf` |
| Figure 2 | Deferred-option illustration with \(d=1/3\) | `scripts/manuscript_figure_2.R` | `figures/generated/manuscript_figure_2.pdf` |
| Figure 3 | Generalized hinge loss | `scripts/manuscript_figure_3.R` | `figures/generated/manuscript_figure_3.pdf` |
| Figure 4 | Overall simulation performance, Settings \((1,2)\) and \((2,2)\) | `scripts/manuscript_figure_4.R` | `figures/generated/manuscript_figure_4.pdf` |
| Figure 5 | A-ERM, ERM, and BO comparison | `scripts/manuscript_figure_5.R` | `figures/generated/manuscript_figure_5.pdf` |
| Figure 6 | PCPT application decisions across \(d\) | `scripts/manuscript_figure_6.R` | `figures/generated/manuscript_figure_6.pdf` |

## Supplementary materials

The supplementary PDF labels its figures “Figure 1” through “Figure 11.” In
this repository they are named S1 through S11 to distinguish them from the
main-manuscript figures.

| Supplementary figure | Caption/topic | Reproduction script | Output |
| --- | --- | --- | --- |
| Figure S1 | Continuous-outcome simulation | `scripts/supplementary_figure_1.R` | `figures/generated/supplementary_figure_1.pdf` |
| Figure S2 | Overall performance, Settings \((2,1)\) and \((2,2)\) | `scripts/supplementary_figure_2.R` | `figures/generated/supplementary_figure_2.pdf` |
| Figure S3 | Deferral metrics, Settings \((2,1)\) and \((2,2)\) | `scripts/supplementary_figure_3.R` | `figures/generated/supplementary_figure_3.pdf` |
| Figure S4 | Overall performance, Settings \((1,1)\) and \((2,1)\) | `scripts/supplementary_figure_4.R` | `figures/generated/supplementary_figure_4.pdf` |
| Figure S5 | Deferral metrics, Settings \((1,1)\) and \((2,1)\) | `scripts/supplementary_figure_5.R` | `figures/generated/supplementary_figure_5.pdf` |
| Figure S6 | Varying IV strength | `scripts/supplementary_figure_6.R` | `figures/generated/supplementary_figure_6.pdf` |
| Figure S7 | Varying covariate dimension | `scripts/supplementary_figure_7.R` | `figures/generated/supplementary_figure_7.pdf` |
| Figure S8 | Simulation under varying confounding strength (point-identifiable/no-\(U\)-to-\(A\) implementation) | `scripts/supplementary_figure_8.R` | `figures/generated/supplementary_figure_8.pdf` |
| Figure S9 | PCPT post-deferral decisions across \(d\) | `scripts/supplementary_figure_9.R` | `figures/generated/supplementary_figure_9.pdf` |
| Figure S10 | PCPT decisions at 90% adherence | `scripts/supplementary_figure_10.R` | `figures/generated/supplementary_figure_10.pdf` |
| Figure S11 | PCPT decisions at 80% adherence | `scripts/supplementary_figure_11.R` | `figures/generated/supplementary_figure_11.pdf` |

## Application-data restriction

Figures 6 and S9–S11 require authorised PCPT saved fit objects. For fits that
do not embed `x` (including some Figure 6 fits, and S10 and S11), provide
either `PCPT_DATA_PATH` or `PCPT_COVARIATE_RDATA`; S9 needs the object
containing the analytic covariates. They are therefore not regenerated from
the public repository alone. The scripts stop with a precise input-data
message when those inputs are unavailable. The public synthetic-data
demonstration remains available through
`scripts/09_run_synthetic_application.R` and
`scripts/figure5_synthetic_application.R`.
