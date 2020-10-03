using RCall
using JellyMe4
using Test

# 'borrowed' from MixedModels.jl until the next major release
using Feather
using Pkg.Artifacts
global TestData = artifact"TestData"
"""
    dataset(nm)

Return the data frame of test data set named `nm`, which can be a `String` or `Symbol`
"""
function dataset(nm::AbstractString)
    path = joinpath(TestData, nm * ".feather")
    if !isfile(path)
        throw(ArgumentError(
            "Dataset \"$nm\" is not available.\nUse MixedModels.datasets() for available names."))
    end
    Feather.read(path)
end
dataset(nm::Symbol) = dataset(string(nm))

"""
    datasets()

Return a vector of names of the available test data sets
"""
datasets() = first.(Base.Filesystem.splitext.(filter(Base.Fix2(endswith, ".feather"), readdir(TestData))))

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
