using RCall
using JellyMe4

using DataFrames
using MixedModels
using Test

const datasets = MixedModels.datasets
dataset(x) = DataFrame(MixedModels.dataset(x))

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


include("merMod.jl")
include("lmerMod.jl")
include("glmerMod.jl")
