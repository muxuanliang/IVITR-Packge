# IV-ITR

### Individualized treatment rules with instrumental variables and a deferred option

Code and reproducibility materials for the accompanying manuscript: R
implementations of the proposed methods, archived simulation results,
figure-generation scripts, and a public synthetic application dataset.

[Manuscript](manuscript.pdf) · [Supplementary materials](sm.pdf) ·
[Figure guide](docs/FIGURE_MAP.md) · [Reproduction instructions](docs/REPRODUCIBILITY.md)

## Start here

| I want to… | Where to start |
| :--- | :--- |
| Recreate the manuscript and supplementary figures | [Figure-by-figure guide](docs/FIGURE_MAP.md) |
| Run the code and display the public figures together | [Executable R Markdown report](reports/reproduce_figures.Rmd) |
| Understand the data-generating models | [Simulation models](docs/SIMULATION_MODELS.md) |
| Rerun simulations or the application demonstration | [Detailed reproduction instructions](docs/REPRODUCIBILITY.md) |
| Check which data are publicly available | [Data availability](docs/DATA_AVAILABILITY.md) |

## Quick start

Use **R 4.4 or later** and run the commands below from the repository root
(the directory containing this README). Archived simulation outputs are
included, so recreating the public figures does not require rerunning the
Monte Carlo experiments.

**1. Install the dependencies and local R package.**

```sh
Rscript scripts/00_install_dependencies.R
```

**2. Generate the public figures and an HTML report.**

```sh
Rscript scripts/render_reproduction_report.R
```

This runs the scripts for manuscript Figures **1–5** and supplementary Figures
**S1–S8**, saves their PDFs under `figures/generated/`, and displays them in
`reproduce_figures.html`. Open that HTML file locally to view the report.
Rendering requires **Pandoc**, available through RStudio or a separate
installation. Keep the HTML file alongside the `figures/` directory because
the report links to the generated PDFs.

## Reproduce an individual figure

Each figure has a script whose name matches its number in the manuscript or
supplement. For example:

```sh
# Manuscript Figure 4: main simulation results
Rscript scripts/manuscript_figure_4.R

# Supplementary Figure 6: simulations under varying IV strength
Rscript scripts/supplementary_figure_6.R
```

These commands produce `figures/generated/manuscript_figure_4.pdf` and
`figures/generated/supplementary_figure_6.pdf`, respectively. Individual figure
scripts do not require Pandoc. The [figure guide](docs/FIGURE_MAP.md) lists all
scripts, outputs, settings, and application-data requirements. The supplement
numbers its figures 1–11; this repository uses S1–S11 to distinguish them from
the main manuscript figures.

## Public application demonstration

A simulated PCPT-like dataset is included with the same variable names and
coding used by the application pipeline. To regenerate the data, fit the
models, and create the plots:

```sh
Rscript scripts/01_prepare_synthetic_application_data.R
Rscript scripts/09_run_synthetic_application.R
Rscript scripts/figure5_synthetic_application.R
```

The fitting step defaults to a 300-subject fit for a smaller demonstration.
See the [data dictionary](data/synthetic/README.md) and
[reproduction instructions](docs/REPRODUCIBILITY.md) for details and options.

> **PCPT application figures:** Manuscript Figure 6 and supplementary Figures
> S9–S11 require authorised PCPT inputs and saved fit objects. The synthetic
> demonstration contains no trial participants and does not reproduce those
> figures exactly. See [data availability](docs/DATA_AVAILABILITY.md).

## Repository contents

| Directory | Contents |
| :--- | :--- |
| [code/IVITR/](code/IVITR/) | R package implementing the methods and comparators |
| [code/simulations/](code/simulations/) | Simulation data generation and estimation |
| [code/figures/](code/figures/) | Plotting and figure assembly |
| [code/application/](code/application/) | Application fitting and plotting code |
| [scripts/](scripts/) | Entry points for figures, simulations, and the application |
| [results/](results/) | Archived Monte Carlo outputs |
| [reports/](reports/) | Executable report source |
| [data/synthetic/](data/synthetic/) | Public synthetic dataset and data dictionary |
| [docs/](docs/) | Models, reproduction instructions, and figure mapping |

Full simulation reruns are computationally intensive and intended for
parallel or HPC execution. The [reproduction instructions](docs/REPRODUCIBILITY.md)
describe task selection and reduced-cost runs. New simulation outputs are
written under `generated/`.

## License

The code is released under the [MIT License](LICENSE).

<details>
<summary>For maintainers: preparing a public release</summary>

Follow the [release checklist](docs/RELEASE_CHECKLIST.md) before publication.
Archive a tagged release in a stable repository such as Zenodo, then add the
repository URL and version-specific DOI to the manuscript's code availability
statement.

</details>
