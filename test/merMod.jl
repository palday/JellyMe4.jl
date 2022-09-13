using RCall, MixedModels, Test
using StatsBase: zscore
const LMM = LinearMixedModel
const GLMM = GeneralizedLinearMixedModel

@testset "merMod" begin
    # this is available in MixedModels.dataset(:sleepstudy) but with different
    # capitalization than in R
    sleepstudy = rcopy(R"sleepstudy")
    jlmm = fit!(LMM(@formula(Reaction ~ 1 + round(Days) + (1 | Subject)), sleepstudy);
                REML=false)
    @testset "bare model" begin
        @test_throws ArgumentError (@rput jlmm)
    end
    @testset "reversed tuple" begin
        jm = (sleepstudy, jlmm)
        @test_throws ArgumentError (@rput jm)
    end
end
