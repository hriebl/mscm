Simulation Study, Scenario 3, Evaluation
================
Hannes Riebl

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(here)
```

    here() starts at /home/hannes/ownCloud/Research/mscm

``` r
library(purrr)

set.seed(1337)

theme_set(theme_minimal())

theme_update(
  text = element_text(family = "TeX Gyre Heros"),
  axis.text = element_text(color = "black", size = rel(1.0)),
  axis.text.x = element_text(margin = margin(2.2, 0, 4.4, 0)),
  strip.text = element_text(color = "black", size = rel(1.0)),
  strip.text.x = element_text(vjust = 0.1),
)

red <- palette()[2]

vars <- c(
  "eta_p0_beta",
  "eta_np0_beta",
  "eta_np0_tau2",
  "gamma",
  "mu_transformed",
  "z",
  "psi",
  "mu"
)

jobs <- list.dirs(here("paper", "simulation-study", "scenario-3", "dont-sync"),
                  recursive = FALSE)

df <- list_rbind(map(jobs, function(job) {
  results <- readRDS(file.path(job, "results.rds"))

  stats <- map(
    .x = vars,
    .f = function(var) {
      truth <- if (var == "mu_transformed") {
        as.vector(log(results$truth[["mu"]]))
      } else {
        as.vector(results$truth[[var]])
      }

      pmean <- as.vector(results$summary$mean[[var]])
      pstddev <- as.vector(results$summary$sd[[var]])

      q <- results$summary$q[[var]]
      commas <- paste0(rep(",", length(dim(q)) - 1), collapse = "")
      q05 <- as.vector(eval(parse(text = paste0("q[1", commas, "]"))))
      q95 <- as.vector(eval(parse(text = paste0("q[3", commas, "]"))))
      covered <- q05 <= truth & truth <= q95

      ess <- as.vector(results$summary$ess_bulk[[var]])
      data.frame(var, truth, pmean, pstddev, covered, ess)
    }
  )

  cbind(as.data.frame(results$param), list_rbind(stats))
}))

df <- df %>%
  mutate(
    nboth = factor(paste0("Sites: ", nsites, "\nSpecies: ", nspecies)),
    var = factor(
      x = var,
      levels = c(
        "eta_p0_beta",
        "eta_np0_beta",
        "eta_np0_tau2",
        "gamma",
        "mu_transformed",
        "z",
        "psi",
        "mu"
      ),
      labels = c(
        "P-β",
        "NP-β",
        "NP-τ²",
        "γ",
        "log(μ)",
        "z",
        "ψ",
        "μ"
      )
    ),
    bias = pmean - truth,
    rmse = sqrt(bias^2)
  )

p <- df %>%
  filter(var == "P-β" | var == "NP-β" | var == "NP-τ²" | var == "γ" | var == "log(μ)") %>%
  group_by(var, nboth) %>%
  slice_sample(n = 5000) %>%
  ungroup()

ggplot(p) +
  geom_jitter(aes(nboth, bias, color = as.factor(nspecies)), alpha = 0.1) +
  geom_boxplot(aes(nboth, bias), fill = NA, outlier.shape = NA) +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(var)) +
  labs(x = "Sample size", y = "Bias") +
  scale_color_discrete(guide = "none") +
  coord_cartesian(ylim = c(-4, 4))
```

![](evaluation_files/figure-commonmark/bias-untransformed-1.png)

``` r
ggplot(p) +
  geom_jitter(aes(nboth, rmse, color = as.factor(nspecies)), alpha = 0.1) +
  geom_boxplot(aes(nboth, rmse), fill = NA, outlier.shape = NA) +
  facet_wrap(vars(var)) +
  labs(x = "Sample size", y = "Root mean squared error (RMSE)") +
  scale_color_discrete(guide = "none") +
  coord_cartesian(ylim = c(0, 4))
```

![](evaluation_files/figure-commonmark/rmse-untransformed-1.png)

``` r
ggplot(p %>% filter(var != "NP-τ²")) +
  geom_point(aes(truth, pmean, color = as.factor(nspecies)), alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(vars(var, nboth), nrow = 1, scales = "free") +
  labs(x = "Truth", y = "Posterior mean") +
  scale_color_discrete(guide = "none")
```

![](evaluation_files/figure-commonmark/truth-untransformed-1.png)

``` r
ggplot(p) +
  geom_jitter(aes(nboth, ess, color = as.factor(nspecies)), alpha = 0.1) +
  geom_boxplot(aes(nboth, ess), fill = NA, outlier.shape = NA) +
  facet_wrap(vars(var)) +
  labs(x = "Sample size", y = "Effective sample size (ESS)") +
  scale_color_discrete(guide = "none")
```

![](evaluation_files/figure-commonmark/ess-untransformed-1.png)

``` r
p <- df %>%
  filter(var == "ψ" | var == "μ") %>%
  group_by(var, nboth) %>%
  slice_sample(n = 5000) %>%
  ungroup()

ggplot(p) +
  geom_jitter(aes(nboth, bias, color = as.factor(nspecies)), alpha = 0.1) +
  geom_boxplot(aes(nboth, bias), fill = NA, outlier.shape = NA) +
  geom_hline(yintercept = 0) +
  facet_wrap(vars(var)) +
  labs(x = "Sample size", y = "Bias") +
  scale_color_discrete(guide = "none") +
  coord_cartesian(ylim = c(-3, 3))
```

![](evaluation_files/figure-commonmark/bias-transformed-1.png)

``` r
ggplot(p) +
  geom_jitter(aes(nboth, rmse, color = as.factor(nspecies)), alpha = 0.1) +
  geom_boxplot(aes(nboth, rmse), fill = NA, outlier.shape = NA) +
  facet_wrap(vars(var)) +
  labs(x = "Sample size", y = "Root mean squared error (RMSE)") +
  scale_color_discrete(guide = "none") +
  coord_cartesian(ylim = c(0, 3))
```

![](evaluation_files/figure-commonmark/rmse-transformed-1.png)

``` r
ggplot(p) +
  geom_point(aes(truth, pmean, color = as.factor(nspecies)), alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(vars(var, nboth), nrow = 1, scales = "free") +
  labs(x = "Truth", y = "Posterior mean") +
  scale_color_discrete(guide = "none")
```

![](evaluation_files/figure-commonmark/truth-transformed-1.png)

``` r
ggplot(p) +
  geom_jitter(aes(nboth, ess, color = as.factor(nspecies)), alpha = 0.1) +
  geom_boxplot(aes(nboth, ess), fill = NA, outlier.shape = NA) +
  facet_wrap(vars(var)) +
  labs(x = "Sample size", y = "Effective sample size (ESS)") +
  scale_color_discrete(guide = "none")
```

![](evaluation_files/figure-commonmark/ess-transformed-1.png)

``` r
p <- df %>%
  filter(var == "z") %>%
  group_by(var, nboth) %>%
  slice_sample(n = 5000) %>%
  ungroup()

ggplot(p) +
  geom_jitter(aes(as.factor(truth), pmean, color = as.factor(nspecies)), alpha = 0.1) +
  facet_wrap(vars(nboth)) +
  labs(x = "Truth", y = "Posterior mean") +
  scale_color_discrete(guide = "none")
```

![](evaluation_files/figure-commonmark/truth-z-1.png)

``` r
df %>%
  group_by(var) %>%
  summarize(bias = mean(bias), .groups = "drop") %>%
  knitr::kable(digits = 3)
```

| var    |   bias |
|:-------|-------:|
| P-β    | -0.001 |
| NP-β   |  0.005 |
| NP-τ²  |  0.312 |
| γ      |  0.000 |
| log(μ) | -0.006 |
| z      |  0.000 |
| ψ      |  0.000 |
| μ      | -0.005 |

``` r
df %>%
  group_by(var) %>%
  summarize(coverage = mean(covered), .groups = "drop") %>%
  knitr::kable(digits = 3)
```

| var    | coverage |
|:-------|---------:|
| P-β    |    0.897 |
| NP-β   |    0.901 |
| NP-τ²  |    0.873 |
| γ      |    0.895 |
| log(μ) |    0.890 |
| z      |    1.000 |
| ψ      |    0.896 |
| μ      |    0.890 |

``` r
df %>%
  group_by(var) %>%
  summarize(
    q01 = quantile(truth, 0.01),
    q99 = quantile(truth, 0.99),
    .groups = "drop"
  ) %>%
  knitr::kable(digits = 3)
```

| var    |    q01 |    q99 |
|:-------|-------:|-------:|
| P-β    | -1.152 |  1.174 |
| NP-β   | -2.312 |  2.327 |
| NP-τ²  |  1.000 |  1.000 |
| γ      | -1.164 |  1.204 |
| log(μ) |  1.432 |  2.933 |
| z      |  0.000 |  1.000 |
| ψ      |  0.120 |  0.883 |
| μ      |  4.185 | 18.781 |
