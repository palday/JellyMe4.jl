module JellyMe4

using RCall
using StatsModels
using MixedModels

export rcopy
export rcopytype
export sexp
export sexpclass

# using Pkg.Artifacts

# function __init__()
#     global TestData = artifact"TestData"
# end

function __init__()
    # should we assume the user is smart enough?
    # reval("library(lme4)")
    if !rcopy(R""" "afex" %in% installed.packages() """)
        error("lme4 must be installed, otherwise what's the point?")
    end

    global LMER = get(ENV, "LMER", "lme4::lmer")
    global AFEX_INSTALLED = rcopy(R""" "afex" %in% installed.packages() """)
    
    nothing
end

function _set_lmer(s)
    global LMER = s
end

function _set_afex_installed(s)
    global AFEX_INSTALLED = s
end


include("utilities.jl")
include("formula.jl")
include("merMod.jl")
include("lmerMod.jl")
include("glmerMod.jl")

end # module
