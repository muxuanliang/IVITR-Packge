#!/usr/bin/env Rscript

# Creates the public PCPT-like synthetic dataset used for the reproducibility
# demonstration. It deliberately contains no participant-level trial data.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else "data/synthetic/pcpt_like.csv"
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

set.seed(20260821)
n <- 5000L
age <- pmin(pmax(round(rnorm(n, 63, 6)), 45), 80)
famhx <- factor(rbinom(n, 1, 0.18), levels = c(0, 1), labels = c("No", "Yes"))
psa_bl <- round(exp(rnorm(n, log(1.1), 0.55)), 2)
smoke <- factor(sample(c("Never", "Former", "Current"), n, replace = TRUE,
                       prob = c(0.42, 0.47, 0.11)))
bmi <- round(pmin(pmax(rnorm(n, 27.2, 4.1), 17), 47), 1)
trt <- ifelse(rbinom(n, 1, 0.5) == 1, "FIN", "PLA")
adherent <- rbinom(n, 1, plogis(1.0 - 0.035 * (age - 63) +
                                  0.20 * (famhx == "Yes") - 0.15 * (bmi - 27) / 4))
adh <- ifelse(trt == "FIN" & adherent == 1, "Adh",
              ifelse(trt == "FIN", "nonAdh", "nonAdh-pill"))
received <- ifelse(trt == "FIN" & adh == "Adh", 1, 0)
linpred <- -1.35 + 0.055 * (age - 63) + 0.50 * log1p(psa_bl) +
  0.30 * (famhx == "Yes") + 0.18 * (smoke == "Current") -
  0.22 * received + 0.09 * received * (psa_bl > 1.4)
high <- rbinom(n, 1, plogis(linpred))
endpoint <- ifelse(high == 1, "High", ifelse(runif(n) < 0.24, "Low", "Healthy"))

synthetic_pcpt <- data.frame(
  age = age, famhx = famhx, psa_bl = psa_bl, smoke = smoke, bmi = bmi,
  endPoint = endpoint, trt = trt, adh = adh
)
utils::write.csv(synthetic_pcpt, output, row.names = FALSE)
message("Wrote ", n, " synthetic records to ", normalizePath(output))
