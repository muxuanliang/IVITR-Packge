# Data availability

## Simulations

No external data are used for the simulation studies. Their complete,
replication-ready Monte Carlo outputs are public in `results/`, and the scripts
that generate every replication are public in `code/simulations/`.

## Application

The manuscript application uses an analytic dataset derived from the Prostate
Cancer Prevention Trial (PCPT). The original participant-level analytic file is
not included in this repository because it is not a file the authors may
redistribute publicly. Prospective users should obtain the source data and any
necessary permissions directly from the data owner, then construct the
analytic file according to their approved data-use terms.

To make the public application demonstration executable without restricted
data, we provide `data/synthetic/pcpt_like.csv`. It is generated with a fixed seed by
`scripts/01_prepare_synthetic_application_data.R`, mimics the variable names,
data types, randomisation, adherence pattern, and binary clinical endpoint used
by the application code, and contains no actual trial participants. It is not
a substitute for the saved PCPT fit objects required to regenerate Figures 6
and S9–S11 exactly.

Before public release, the authors should confirm this statement with the PCPT
data owner and replace it if the analytic file can instead be deposited in a
public, stable repository.
