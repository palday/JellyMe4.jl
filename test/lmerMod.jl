using RCall, MixedModels, Test
using StatsBase: zscore
using Tables: columntable
const LMM = LinearMixedModel
const GLMM = GeneralizedLinearMixedModel

@testset "lmerMod" begin
    reval("""
    if(!require(lme4)){
        # use tmp for tests if lme4 isn't available
        .libPaths("/tmp")
        lib <- .libPaths()[1L]
        install.packages("lme4",repos="https://cloud.r-project.org", libs=lib)
        library(lme4)
    }""")

    # this is available in MixedModels.dataset(:sleepstudy) but with different
    # capitalization than in R
    sleepstudy = rcopy(R"sleepstudy")

    @testset "get lmerMod" begin
        ### from R ###
        jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + (1 + Days|Subject)),sleepstudy), REML=false)
        rlmm = rcopy(R"m <- lmer(Reaction ~ 1 + Days + (1 + Days|Subject),sleepstudy,REML=FALSE)")

        @test jlmm.θ ≈ rlmm.θ atol=0.001
        @test objective(jlmm) ≈ objective(rlmm) atol=0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001

        jlmm = fit!(jlmm, REML=true)
        rlmm = rcopy(R"update(m, REML=TRUE)")

        @test jlmm.θ ≈ rlmm.θ atol=0.001
        @test objective(jlmm) ≈ objective(rlmm) atol=0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001

        # single scalar RE is different because RCall converts the single-entry θ
        # to a scalar value, which we have to correct
        # of course, this isn't a problem in the other direction, because
        # the scalar-vector distinction for θ is missing in R
        @testset "scalar RE" begin
            jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + (1|Subject)),sleepstudy), REML=false)
            rlmm = rcopy(R"m <- lmer(Reaction ~ 1 + Days + (1|Subject),sleepstudy,REML=FALSE)")

            @test jlmm.θ ≈ rlmm.θ atol=0.001
            @test objective(jlmm) ≈ objective(rlmm) atol=0.001
            @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001
        end
    end

    @testset "put lmerMod" begin
        ### from Julia ###
        jlmm = LMM(@formula(Reaction ~ 1 + Days + (1 + Days|Subject)),sleepstudy)
        jm = Tuple([jlmm, sleepstudy])
        # unfitted model
        @test_throws ArgumentError @rput jm
        fit!(jlmm, REML=true)
        @rput jm
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
        @test rcopy(R"REMLcrit(jm)") ≈ objective(jlmm)

        fit!(jlmm, REML=false)
        jm = Tuple([jlmm, sleepstudy])
        @rput jm
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
        @test rcopy(R"deviance(jm)") ≈ objective(jlmm)

        @testset "columntable" begin
            jm = Tuple([jlmm, columntable(sleepstudy)]);
            @rput jm;
        end

        @testset "transformations" begin
            sleepstudy[!,:Days2] = sleepstudy.Days .+ 1
            @rput sleepstudy;
            R"m <- lmer(log10(Reaction) ~ 1 + log(Days2) + (1 + log(Days2)|Subject),sleepstudy,REML=FALSE)"
            jlmm = fit!(LMM(@formula(log10(Reaction) ~ 1 + log(Days2) + (1 + log(Days2)|Subject)),sleepstudy), REML=false)
            jm = Tuple([jlmm, sleepstudy])
            @rput jm;
            @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
            @test rcopy(R"deviance(jm)") ≈ objective(jlmm)

            jlmm = fit!(LMM(@formula(Reaction ~ 1 + round(Days) + (1|Subject)),sleepstudy), REML=false)
            jm = Tuple([jlmm, sleepstudy]);

            @test_throws ArgumentError (@rput jm)
        end

    end
end
