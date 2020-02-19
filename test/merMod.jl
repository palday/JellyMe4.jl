using RCall, MixedModels, Test
using StatsBase: zscore
const LMM = LinearMixedModel
const GLMM = GeneralizedLinearMixedModel

@testset "merMod" begin
    reval("""
    if(!require(lme4)){
        # use tmp for tests if lme4 isn't available
        .libPaths("/tmp")
        lib <- .libPaths()[1L]
        warning("installing lme4")
        warning(lib)
        install.packages("lme4",repos="https://cloud.r-project.org", libs=lib)
        library(lme4)
    }""")

    # this is available in MixedModels.dataset(:sleepstudy) but with different
    # capitalization than in R
    sleepstudy = rcopy(R"sleepstudy")
    jlmm = fit!(LMM(@formula(Reaction ~ 1 + round(Days) + (1|Subject)),sleepstudy), REML=false)
    @testset "bare model" begin
        @test_throws ArgumentError (@rput jlmm)
    end
end
