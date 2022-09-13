module JellyMe4

using RCall
using StatsModels
using RegressionFormulae
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
    if !rcopy(R""" "lme4" %in% installed.packages() """)
        @error "lme4 not found"
    end

    global LMER = get(ENV, "LMER", "lme4::lmer")
    global AFEX_INSTALLED = rcopy(R""" "afex" %in% installed.packages() """)

    return nothing
end

function _set_lmer(s)
    return global LMER = s
end

function _set_afex_installed(s)
    return global AFEX_INSTALLED = s
end

const MERCONTROL_OPTIONS = ["""optimizer="nloptwrap" """,
                            """calc.derivs=FALSE""",
                            """check.nobs.vs.rankZ = "warning" """,
                            """check.nobs.vs.nlev = "warning" """,
                            """check.nlev.gtreq.5 = "ignore" """,
                            """check.nlev.gtr.1 = "warning" """,
                            """check.nobs.vs.nRE= "warning" """,
                            """check.rankX = "message+drop.cols" """,
                            """check.scaleX = "warning" """,
                            """check.formula.LHS = "stop" """]

include("utilities.jl")
include("formula.jl")
include("merMod.jl")
include("lmerMod.jl")
include("glmerMod.jl")

end # module
