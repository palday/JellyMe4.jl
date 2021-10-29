# JellyMe4

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://palday.github.io/JellyMe4.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://palday.github.io/JellyMe4.jl/dev)
[![Build Status](https://travis-ci.com/palday/JellyMe4.jl.svg?branch=master)](https://travis-ci.com/palday/JellyMe4.jl)
-->
[![Test Status](https://github.com/palday/JellyMe4.jl/workflows/CI/badge.svg)](https://github.com/palday/JellyMe4.jl/actions)
[![Codecov](https://codecov.io/gh/palday/JellyMe4.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/palday/JellyMe4.jl)
[![DOI](https://zenodo.org/badge/240243868.svg)](https://zenodo.org/badge/latestdoi/240243868)

<!-- [![Build Status](https://ci.appveyor.com/api/projects/status/github/palday/JellyMe4.jl?svg=true)](https://ci.appveyor.com/project/palday/JellyMe4-jl) -->

* [Purpose](#purpose)
* [Installation](#installation)
* [Basic Usage](#basic-usage)
  + [Fitting a model in R and moving it to Julia](#fitting-a-model-in-r-and-moving-it-to-julia)
  + [Fitting a model in Julia and moving it to R](#fitting-a-model-in-julia-and-moving-it-to-r)
* [Limitations and warnings](#limitations-and-warnings)
* [Alternative `lmer`](#alternative-lmer)
  + [Julia to R](#julia-to-r)
  + [R to Julia](#r-to-julia)
* [Where does the name come from?](#where-does-the-name-come-from)
* [Acknowledgements](#acknowledgements)


## Purpose
One of the difficulties in transitioning to a new programming language is not just learning how to do things in the new language, but the difference in the package ecosystem. `RCall` helps with both aspects when moving from R to Julia, at least when it comes to basic data manipulation. `JellyMe4` takes advantage of `RCall`'s extensibility to provide a way to transfer mixed-effects models fit between R's `lme4` and Julia's `MixedModels`. This means that it is now possible to fit a model in Julia, but then take advantage of existing R packages for examining the model such as `car` and `effects`.

## Installation

`JellyMe4` is now registered in the Julia package registry and can be installed by
```julia
julia> using Pkg
julia> Pkg.add("JellyMe4")

```
or, in the pkg REPL,
```julia
(@v1.4) pkg> add JellyMe4
```

To get the latest pre-release features, you can install the development version:
```julia
(@v1.4) pkg> add JellyMe4#master
```

Generally speaking, the development version should work, but especially until version 1.0, there is no guarantee that there won't be breaking changes compared to the latest release.

## Basic Usage

### Fitting a model in R and moving it to Julia

```julia
julia> using MixedModels, RCall, JellyMe4
julia> R"library(lme4)"
┌ Warning: RCall.jl: Loading required package: Matrix
└ @ RCall ~/.julia/packages/RCall/g7dhB/src/io.jl:113
RObject{StrSxp}
[1] "lme4"      "Matrix"    "stats"     "graphics"  "grDevices" "utils"
[7] "datasets"  "methods"   "base"

julia> R"summary(m <- lmer(Reaction ~ 1 + Days + (1+Days|Subject),sleepstudy,REML=FALSE))"
RObject{VecSxp}
Linear mixed model fit by maximum likelihood  ['lmerMod']
Formula: Reaction ~ 1 + Days + (1 + Days | Subject)
   Data: sleepstudy

     AIC      BIC   logLik deviance df.resid
  1763.9   1783.1   -876.0   1751.9      174

Scaled residuals:
    Min      1Q  Median      3Q     Max
-3.9416 -0.4656  0.0289  0.4636  5.1793

Random effects:
 Groups   Name        Variance Std.Dev. Corr
 Subject  (Intercept) 565.48   23.780
          Days         32.68    5.717   0.08
 Residual             654.95   25.592
Number of obs: 180, groups:  Subject, 18

Fixed effects:
            Estimate Std. Error t value
(Intercept)  251.405      6.632  37.907
Days          10.467      1.502   6.968

Correlation of Fixed Effects:
     (Intr)
Days -0.138


julia> @rget m
Linear mixed model fit by maximum likelihood
 Reaction ~ 1 + Days + (1 + Days | Subject)
   logLik   -2 logLik     AIC        BIC
 -875.96967 1751.93934 1763.93934 1783.09709

Variance components:
            Column    Variance   Std.Dev.    Corr.
Subject  (Intercept)  565.476967 23.7797596
         Days          32.681785  5.7167985  0.08
Residual              654.945706 25.5919070
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
             Estimate  Std.Error   z value  P(>|z|)
───────────────────────────────────────────────────
(Intercept)  251.405     6.63212  37.9072    <1e-99
Days          10.4673    1.50223   6.96783   <1e-11
───────────────────────────────────────────────────
julia> typeof(m)
LinearMixedModel{Float64}
```

### Fitting a model in Julia and moving it to R
```julia
julia> using MixedModels, RCall, JellyMe4
julia> machines = rcopy(R"nlme::Machines")
54×3 DataFrames.DataFrame
│ Row │ Worker       │ Machine      │ score   │
│     │ Categorical… │ Categorical… │ Float64 │
├─────┼──────────────┼──────────────┼─────────┤
│ 1   │ 1            │ A            │ 52.0    │
│ 2   │ 1            │ A            │ 52.8    │
│ 3   │ 1            │ A            │ 53.1    │
│ 4   │ 2            │ A            │ 51.8    │
│ 5   │ 2            │ A            │ 52.8    │
│ 6   │ 2            │ A            │ 53.1    │
⋮
│ 48  │ 4            │ C            │ 64.0    │
│ 49  │ 5            │ C            │ 72.1    │
│ 50  │ 5            │ C            │ 72.0    │
│ 51  │ 5            │ C            │ 71.1    │
│ 52  │ 6            │ C            │ 62.0    │
│ 53  │ 6            │ C            │ 61.4    │
│ 54  │ 6            │ C            │ 60.5    │

julia> m = fit(MixedModel, @formula(score ~ 1 + Machine + (1 + Machine|Worker)), machines)
Linear mixed model fit by maximum likelihood
 score ~ 1 + Machine + (1 + Machine | Worker)
   logLik    -2 logLik      AIC         BIC
 -108.208914  216.417828  236.417828  256.307668

Variance components:
            Column     Variance   Std.Dev.   Corr.
Worker   (Intercept)  13.81760447 3.7172039
         Machine: B   28.68496515 5.3558347  0.49
         Machine: C   11.24097825 3.3527568 -0.36  0.30
Residual               0.92463436 0.9615791
 Number of obs: 54; levels of grouping factors: 6

  Fixed-effects parameters:
───────────────────────────────────────────────────
             Estimate  Std.Error   z value  P(>|z|)
───────────────────────────────────────────────────
(Intercept)  52.3556     1.53437  34.1218    <1e-99
Machine: B    7.96667    2.20988   3.60502   0.0003
Machine: C   13.9167     1.40579   9.89956   <1e-22
───────────────────────────────────────────────────
julia> # LinearMixedModel doesn't keep of the original dataframe,
julia> # so we need to package it up
julia> m_machines = (m, machines);
julia> @rput m_machines;
julia> R"summary(m_machines)"
RObject{VecSxp}
Linear mixed model fit by maximum likelihood  ['lmerMod']

     AIC      BIC   logLik deviance df.resid
   236.4    256.3   -108.2    216.4       44

Scaled residuals:
     Min       1Q   Median       3Q      Max
-2.40773 -0.51890  0.03227  0.45598  2.54091

Random effects:
 Groups   Name        Variance Std.Dev. Corr
 Worker   (Intercept) 13.8176  3.7172
          MachineB    28.6850  5.3558    0.49
          MachineC    11.2410  3.3528   -0.36  0.30
 Residual              0.9246  0.9616
Number of obs: 54, groups:  Worker, 6

Fixed effects:
            Estimate Std. Error t value
(Intercept)   52.356      1.534  34.122
MachineB       7.967      2.210   3.605
MachineC      13.917      1.406   9.900

Correlation of Fixed Effects:
         (Intr) MachnB
MachineB  0.463
MachineC -0.374  0.301
```

Note that when moving models from Julia to R, you may see warnings like:
```julia
┌ Warning: RCall.jl: Warning in optwrap(optimizer, devfun, start, rho$lower, control = control,  :
│   convergence code 5 from nloptwrap: NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached.
```

This is expected because we don't allow optimization to stop naturally when moving models back forth, but instead force optimization to stop *because we already found an optimum in the other language*.

## Limitations and warnings

This is alpha software. It has some functionality that should work well for common use cases and even a testsuite, but this testsuite depends on two different software environments (R and Julia) and can afoul of all sorts of nasty version interactions. The testuite only tests against a single version of R and the current version of lme4. In other words, even for the parts that do work for me, they may not work for you.

Parts of the API aren't as nice as I would like them to be (especially in the Julia to R direction) and may change if I figure out something nicer.
This package generally follows the Julia community standard [SemVer](https://semver.org/).
As this package has not yet hit 1.0, we follow these general conventions:
  - Non-exported functions may change or disappear even with a patch release. Exported functions should otherwise only increase in functionality between patch versions; supporting additional syntax is not considered a breaking change. In other words, errors and exceptions disappearing as more conversions are supported is not breaking. Small improvements will be accompanied by a patch version bump.
  - Large improvements and potentially breaking changes will be accompanied by a minor version bump.

Only a subset of all the options available in `MixedModels` and `lme4` are supported. Unsupported things should break in an obvious way, but here's a list of things that are more commonly used but may break non obviously:
- ~~**custom contrast coding**. If you really need this and you know what you're doing, then you can set up your own numeric variables representing appropriate contrasts, as numeric variables survive the transition without difficulty.~~ This should work if you're not doing anything weird with your interaction terms.
- ~~**advanced models in either language** (e.g., `zerocorr!`, `fulldummy` in Julia or `||` in R, which are not completely synonymous anyway).~~ This should largely work, but there are a few edge cases, especially when combining these, that might not work. `dummy` in R will not be translated -- create your indicator variables ahead of time and not as part of the formula.  `zerocorr` will work with continuous variables and categorical variables with two levels using `lme4` on the R side. Categorical variables with more than two levels require the R package [`afex`](https://cran.r-project.org/web/packages/afex) to be installed for `zerocorr` to be translated correctly. See [Alternative lmer](#alternative--lmer-) below. Using `||` in R and translating to Julia will currently fail with categorical variables.
- **fancy transformations within model formulae**, especially those that have different names (e.g. `scale()` in R). If in doubt, pre-transform your data in the dataframe before fitting. A few transformations are supported (e.g. `log`, `exp`) and you should receive an error for unsupported transformations, but it's always a good idea to compare the estimates of the original and copied models.
- **missing data** is handled differently in R and Julia, both in its representation and when it's eliminated. If in doubt, remove missing data before fitting models for consistency.
- ~~**interactions specified in R with `^`**. This part of the parsing is handled by RCall and winds up getting translated into simple exponentiation instead of "interactions below a certain order".  Consider using the [`simplify.formula`](https://www.rdocumentation.org/packages/MuMIn/versions/1.9.5/topics/Formula%20manipulation) function from the [`MuMIn` package](https://cran.r-project.org/web/packages/MuMIn/index.html) in R to simplify your formula to primitive operations before calling lme4.~~ This should now work, by using RegressionFormulae.jl.
- **other R-specific extensions to Wilkinson-Roger notation**. This includes `%in%` for nesting, (`/` for nesting, with reversed order should work because lme4 and MixedModels both handle that specially), `I()` for literal interpretation of arithmetic, `offset()`.
- **rank deficiency in the fixed effects** is handled differently in R and Julia. If in doubt, [remove extraneous columns/predictors](https://juliastats.org/MixedModels.jl/dev/rankdeficiency/#Rank-deficiency-in-mixed-effects-models) before fitting models for consistency.
- **GLMMs with dispersion parameters**. This should error when unsupported model types are encountered. This is intentionally not supported at the moment because there are [known problems with these models in MixedModels.jl](JuliaStats/MixedModels#291).
- **getting binomial GLMMs from R with any of the alternative specifications**, e.g. `cbind(successes, failures)` or `successes/total`. You will receive a gentle reminder in the form of an error. Note that even once this is supported, you're subject to the constraints above.
- **getting GLMMs with `quasi` families from R**. The corresponding distributions aren't in Distributions.jl and are only "quasi" in the GLM/GLMM framework.
- **a number of scratch variables in R prefixed by `jellyme4_` are created.** We need scratch variables when moving things to R, so we use  identifiers beginning with `jellyme4_`.

Finally, to work its magic, this package hooks into `MixedModels` and `lme4` internals. `lme4`'s internals are quite stable now and this is taken advantage of by several other R packages (`ordinal`, `gamm4`, `robustlmm`). `MixedModels` are also fairly stable for the core things, but not everything. I am a contributor to `MixedModels`, so I generally know when something is going to change, but things may still break with changing `MixedModels` versions, especially for major version changes.

## Alternative `lmer`

### Julia to R
By default, JellyMe4 uses `[g]lmer` from `lme4`, as this is the closest equivalent to `MixedModels` (and was in no small part written by the same individual). It is, however, possible to set a different `lmer` implementation. Notably, the `afex` and `lmerTest` packages extend `lmer` in ways that are often appealing. If you set the Julia environment variable `ENV["LMER"]` before loading JellyMe4, then JellyMe4 will use the alternative specification. For example:
```julia
# note the use of R's package::function syntax
julia> ENV["LMER"] = "lmerTest::lmer" # if you really need Satterthwaite or Kenward-Roger ddf
julia> ENV["LMER"] = "afex::lmer_alt" # for better handling of the || syntax
julia> using RCall, JellyMe4
```
There is currently no public API for changing this after loading the package. For `lmerTest`, it is possible to modify the resultant `lme4::lmerMod` object to be an `lmerTest::merModLmerTest` after the fact:
It is also possible to do this conversion after the fact :
```R
R> mod = as(mod, "merModLmerTest")
```
(Note that this conversion is only possible the LMM case. For GLMM, the degrees of freedom are not an issue.)
If `afex` is installed and you attempt to convert a `zerocorr` model with a categorical variable with many levels, then JellyMe4 will swap itself to use `afex::lmer_alt` for the correct handling of these. If `afex` is not installed, then JellyMe4 can still convert two-level factors and continuous variables using the machinery built into `lme4`, but will error for multi-level factors. Once JellyMe4 has swapped to using `afex::lmer_alt`, it will continue doing so for the remainder of the session, but it will only swap when it first encounters a multi-level factor.

### R to Julia

Models fit with the R package [`lmerTest`](https://cran.r-project.org/web/packages/lmerTest/) should work without difficulty as `lmerTest` does not modify the fitting process nor the relevant internal structure of the model of `lme4`.
Note that the functionality related to denominator degrees of freedom (e.g. Satterthwaite or Kenward-Roger approximations) are not currently implemented in `MixedModels` and thus this functionality will be lost.

Models fit with the R package [`afex`](https://cran.r-project.org/web/packages/afex/) should also largely be compatible, as `afex` works by providing a more convenient interface to `lme4`.
Note however that if any of `afex`'s behind the scenes rewrite-rules invoke features not yet supported, then the resulting models will not work with JellyMe4.
In particular, `afex` directly modifies the model matrix to fix the limitation in the `||` syntax and this will almost definitely not work in JellyMe4.
This could probably be made to work, but it seems needlessly complicated given that the assumption is that most models are easier to fit with `MixedModels` than `lme4` and that the package is primarily for moving models from Julia to R.
In any case, please check your model summary on both sides to make sure they line up!

## Where does the name come from?

A really bad pun.

If you take `lme4` and add a J for Julia and try to say the resulting string quickly, it sounds a bit like the package name.

Sorry, not sorry.

## Acknowledgements

The development of this package was supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".
