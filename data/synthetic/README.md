# Public synthetic PCPT-like data

`pcpt_like.csv` is generated deterministically by
`scripts/01_prepare_synthetic_application_data.R` using seed `20260821`.
It is simulated data, not PCPT participant data.

| Variable | Description |
| --- | --- |
| `age` | Baseline age in years. |
| `famhx` | Family-history indicator (`No`/`Yes`). |
| `psa_bl` | Baseline PSA-like continuous measurement. |
| `smoke` | Smoking category. |
| `bmi` | Body-mass index. |
| `endPoint` | `Healthy`, `Low`, or `High` endpoint category. |
| `trt` | Randomised assignment: `FIN` or `PLA`. |
| `adh` | Adherence category used to construct treatment received. |
