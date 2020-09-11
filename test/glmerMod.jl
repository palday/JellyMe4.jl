using RCall, MixedModels, Test
using StatsBase: zscore
using Tables: columntable
const LMM = LinearMixedModel
const GLMM = GeneralizedLinearMixedModel

logistic(x)  = 1 / (1 + exp(-x))

@testset "glmerMod" begin
    ### from Julia ###
    @testset "put glmerMod" begin
        # TODO: remove upon next MixedModels.jl release
        dat = dataset(:verbagg);
        jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)+(1|item)),
                    dat, Bernoulli());

        jm = (jlmm, dat);
        # unfitted model
        @test_throws ArgumentError @rput jm
        fit!(jlmm, fast=false)
        @rput jm;
        # @test_warn Regex(".*categorical.*") @rput jm;
        # note the really high tolerances
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.1
        @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.5

        jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)+(1|item)),
                    dat, Bernoulli());
        fit!(jlmm, fast=true)
        jm = (jlmm, dat)
        @rput jm;
        # note the really high tolerances
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.1
        @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.5;

        @testset "columntable" begin
            jm = (jlmm, columntable(dat));
            @rput jm;
        end

        @testset "contrasts" begin
            jlmm = fit(MixedModel, @formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)+(1|item)),
                        dat, Bernoulli(), contrasts=Dict(:gender => EffectsCoding(),
                                                         :btypes => EffectsCoding()));
            jm = (jlmm, dat);
            @rput jm;
            @test fixef(jlmm) ≈  rcopy(R"fixef(jm)");
        end

    end
end
