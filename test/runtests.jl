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


include("merMod.jl")
include("lmerMod.jl")
include("glmerMod.jl")
