module JellyMe4

using StatsModels

export rcopy
export rcopytype
export sexp
export sexpclass

# using Pkg.Artifacts

# function __init__()
#     global TestData = artifact"TestData"
# end

include("formula.jl")
include("merMod.jl")
include("lmerMod.jl")
include("glmerMod.jl")

end # module
