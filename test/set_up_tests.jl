
using DataFrames
using GLM
using JellyMe4
using MixedModels
using RCall
using Test

using GLM: Link
using JellyMe4: _set_lmer, _set_afex_installed
using MixedModelsDatasets: MixedModelsDatasets, datasets
using StatsBase: zscore
using StatsModels: SeqDiffCoding
using Tables: columntable

const GLMM = GeneralizedLinearMixedModel
const LMM = LinearMixedModel

dataset(x) = DataFrame(MixedModelsDatasets.dataset(x))

# this should only occur on CIs
# the Linux CI installs lme4 via apt beforehand
# and there are binary builds for Mac and Windows
reval("""
if(!require(lme4)){
    # use tmp for tests if lme4 isn't available
    .libPaths("/tmp")
    lib <- .libPaths()[1L]
    warning("installing lme4")
    warning(lib)
    install.packages("lme4",repos="https://cloud.r-project.org", libs=lib, type="binary")
    library(lme4)
}
if(!require(afex)){
    .libPaths("/tmp")
    lib <- .libPaths()[1L]
    install.packages("afex",repos="https://cloud.r-project.org", libs=lib, , type="binary")
}
""")


logistic(x) = 1 / (1 + exp(-x))
