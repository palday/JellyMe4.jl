module JellyMe4

# if RCall is available, then so is DataFrames
# so don't feel bad about this dependency
using DataFrames: DataFrames,
                  DataFrame,
                  transform!

using CategoricalArrays: categorical,
                         CategoricalArray

using Distributions: Distribution,
                     Bernoulli,
                     Binomial,
                     Gamma,
                     Gaussian,
                     InverseGaussian,
                     Poisson
using GLM: CauchitLink,
           CloglogLink,
           IdentityLink,
           InverseLink,
           Link,
           LogitLink,
           LogLink,
           ProbitLink,
           SqrtLink

using LinearAlgebra: pinv

using MixedModels: MixedModels,
                   MixedModel,
                   getθ,
                   GeneralizedLinearMixedModel,
                   LinearMixedModel,
                   nranef,
                   pirls!,
                   setβθ!,
                   setθ!,
                   updateL!

using RCall: RCall,
             @rget,
             @rput,
             @R_str,
             LangSxp,
             protect,
             rcopy,
             rcopytype,
             reval,
             RObject,
             RClass,
             S4Sxp,
             sexp,
             sexpclass,
             unprotect

using RegressionFormulae: RegressionFormulae

using StatsModels: StatsModels,
                   @formula,
                   CategoricalTerm,
                   coefnames,
                   HypothesisCoding,
                   InterceptTerm,
                   InteractionTerm

using Tables: columntable,
              ColumnTable

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

const MERCONTROL_OPTIONS = ["""calc.derivs=FALSE""",
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
