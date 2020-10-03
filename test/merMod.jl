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
        # this should only occur on CIs
        # the Linux CI installs lme4 via apt beforehand
        # and there are binary builds for Mac and Windows
        install.packages("lme4",repos="https://cloud.r-project.org", libs=lib, type="binary")
        library(lme4)
    }
    if(!require(afex)){
        .libPaths("/tmp")
        lib <- .libPaths()[1L]
        install.packages("afex",repos="https://cloud.r-project.org", libs=lib, , type="binary")
    }
    """)

    # this is available in MixedModels.dataset(:sleepstudy) but with different
    # capitalization than in R
    sleepstudy = rcopy(R"sleepstudy")
    jlmm = fit!(LMM(@formula(Reaction ~ 1 + round(Days) + (1|Subject)),sleepstudy), REML=false)
    @testset "bare model" begin
        @test_throws ArgumentError (@rput jlmm)
    end
    @testset "reversed tuple" begin
        jm = (sleepstudy, jlmm);
        @test_throws ArgumentError (@rput jm)
    end
end
