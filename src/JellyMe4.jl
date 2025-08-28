module JellyMe4

using RCall
using StatsModels
using RegressionFormulae
using MixedModels

import MixedModels: GeneralizedLinearMixedModel,
                    pirls!,
                    setβθ!,
                    setθ!
import Distributions: Distribution,
                      Bernoulli,
                      Binomial,
                      Gamma,
                      Gaussian,
                      InverseGaussian,
                      Poisson

import GLM: Link,
            CauchitLink,
            CloglogLink,
            LogitLink,
            IdentityLink,
            InverseLink,
            LogLink,
            ProbitLink,
            SqrtLink

import RCall: rcopy,
              RClass,
              rcopytype,
              reval,
              S4Sxp,
              sexp,
              protect,
              unprotect,
              sexpclass,
              @rput,
              @rget
# if RCall is available, then so is DataFrames
using DataFrames
import CategoricalArrays: CategoricalArray
import Tables: columntable, ColumnTable

import MixedModels: LinearMixedModel,
                    setθ!,
                    updateL!
import RCall: rcopy,
              RClass,
              rcopytype,
              reval,
              S4Sxp,
              sexp,
              protect,
              unprotect,
              sexpclass,
              @rput,
              @rget,
              @R_str

# if RCall is available, then so is DataFrames
import DataFrames: DataFrame
import Tables: ColumnTable

import RCall: rcopy,
              RObject,
              LangSxp,
              @rput,
              @rget,
              @R_str

import MixedModels: MixedModel
import Tables: ColumnTable
import DataFrames: DataFrame
import RCall: rcopy,
              RClass,
              rcopytype,
              reval,
              S4Sxp,
              sexp,
              protect,
              unprotect,
              sexpclass,
              @rput,
              @rget

using CategoricalArrays
import LinearAlgebra: pinv
using DataFrames
using MixedModels: getθ, nranef
import RCall: rcopy,
              reval,
              @R_str

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
