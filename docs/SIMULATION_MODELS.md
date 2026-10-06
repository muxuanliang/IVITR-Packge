---
output:
  pdf_document: default
  html_document: default
---
# Simulation data-generating models

Let \(X=(X_1,\ldots,X_p)\), where the components are independent
\(\mathrm{Uniform}(-1,1)\), let \(U\sim\mathrm{Uniform}(-1,1)\), and let
\(Z^*\in\{-1,1\}\) be randomised with probability one half. The observed
treatment and outcome are binary. The scripts recode binary variables to
\(\{-1,1\}\) before fitting.

## Main simulation

For the main setting, \(p=5\), \(\lambda=\xi=\delta=\eta\), and

\[
 A\mid X,U,Z^*\sim\mathrm{Bernoulli}\{\operatorname{expit}(8Z^*+X_1-7X_2+
 \lambda(1+X_1)U)\}.
\]

Conditional on \((X,U,A)\),
\[
 Y\sim\mathrm{Bernoulli}\{\operatorname{expit}(g_1(X,U)+g_2(X,U)A)\},
\]
where

\[
g_1(X,U)=
\begin{cases}
1-X_1-X_2+\xi U,&\text{setting 1}=1,\\
1-X_1^2+X_2^2+\xi X_1X_2U,&\text{setting 1}=2,
\end{cases}
\]

and

\[
g_2(X,U)=
\begin{cases}
0.442(1-X_1-X_2+\delta U),&\text{setting 2}=1,\\
X_2-0.25\delta X_1^2+4I(X_2>0)-6I(X_2\le0),&\text{setting 2}=2.
\end{cases}
\]

The main grid uses \(n\in\{300,500\}\), \(\eta\in\{0.5,1.5\}\), both
setting indicators in \(\{1,2\}\), and
\(d\in\{0.5,0.4,0.3,0.2,0.1\}\).

## Additional simulations under varying design features

These additional simulations retain the same outcome model and vary one
first-stage feature at a time: IV strength, covariate dimension, or the
strength of unmeasured confounding in treatment assignment. They otherwise
use setting \((2,2)\) with \(\eta=1.5\). The mild setting removes the
\(U\)-to-\(A\) term and uses
\(\operatorname{logit}P(A=1\mid Z^*,X)=2Z^*+X_1-X_2\). See
`code/simulations/test_variant.R` and `code/simulations/test_mild.R` for the
exact coefficients and grids.

## Identifiable setting

The identifiable setting makes the minimal first-stage change of removing the
unmeasured confounder from treatment assignment while retaining the main
manuscript first-stage coefficients:

\[
 A\mid X,Z^*\sim\mathrm{Bernoulli}\{\operatorname{expit}(8Z^*+X_1-7X_2)\}.
\]

The outcome model and the four combinations of settings 1 and 2 are exactly as
above, with \(p=5\) and \(\xi=\delta=\eta\). Consequently the conditional
treatment effect is identified by the Wald-type comparison without deferral.
The proposed method is still evaluated over the same five values of \(d\).

For this setting, the scripts report both the worst-case loss computed from
the population partial-identification bounds and the loss using the true CATE.
The code checks at every replication that the former is no smaller than the
latter under the implemented reparameterisation of \(d\).
