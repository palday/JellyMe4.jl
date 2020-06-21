using RCall, MixedModels, Test
using StatsBase: zscore
using Tables: columntable
const LMM = LinearMixedModel
const GLMM = GeneralizedLinearMixedModel

logistic(x)  = 1 / (1 + exp(-x))

@testset "glmerMod" begin
    ### from Julia ###
    @testset "put glmerMod" begin
        @testset "Bernoulli" begin
            # NOTE: For Bernoulli, some of the tests are doubled in order to
            #       check nAQG and fast options
            # TODO: remove upon next MixedModels.jl release
            dat = dataset(:verbagg);
            # we max out at 1 scalre RE for the nAGQ tests
            jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)),
                        dat, Bernoulli());

            @testset "unfitted model" begin
                jm = Tuple([jlmm, dat]);
                # unfitted model
                @test_throws ArgumentError @rput jm
            end

            @testset "nAGQ" begin
                fit!(jlmm; fast=false, nAGQ=9)
                @rput jm;
                # @test_warn Regex(".*categorical.*") @rput jm;

                @test  rcopy(R"""jm@devcomp$dims["nAGQ"]""") == jm.optstum.nAGQ
                # note the really high tolerances
                @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.1
                @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.1
            end

            @testset "Laplace" begin
                jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)+(1|item)),
                            dat, Bernoulli());
                fit!(jlmm; fast=true)
                jm = Tuple([jlmm, dat])
                @rput jm;

                # note the really high tolerances
                @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.1
                @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.5;
            end
        end

        @testset "columntable" begin
            jm = Tuple([jlmm, columntable(dat)]);
            @rput jm;
        end

        @testset "Binomial" begin
            # TODO: remove upon next MixedModels.jl release
            dat  = dataset(:cbpp);
            dat.rate = dat.incid ./ dat.hsz;
            jlmm = fit(MixedModel, @formula(rate ~ 1 + period + (1|herd)),
                      dat, Binomial(); wts=float(dat.hsz));

            jm = (jlmm, dat);
            fit!(jlmm, fast=true)
            @rput jm;
            # @test_warn Regex(".*categorical.*") @rput jm;
            @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.001
            @test_broken rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.001
            # this is where the weights have their biggest impact
            @test stderr(jlmm) ≈ rcopy(R"""coef(summary(jm))[,"Std. Error"]""") atol=0.001
            @test rcopy(R"logLik(jm)") ≈ loglikelihood(jlmm) atol=0.001
        end

        @testset "Poisson" begin
            # TODO: remove upon next MixedModels.jl release
            center(v::AbstractVector) = v .- (sum(v) / length(v))
            dat = dataset(:grouseticks);
            dat.ch = center(dat.height);
            jlmm = GLMM(@formula(ticks ~ 1+year+ch+ (1|index) + (1|brood) + (1|location)),
                        dat, Poisson());

            jm = (jlmm, dat);
            fit!(jlmm, fast=true) # problems with this one in fast=false
            @rput jm;
            # @test_warn Regex(".*categorical.*") @rput jm;
            # note the really high tolerances
            @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.001
            @test_broken rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.5
            @test rcopy(R"logLik(jm)") ≈ loglikelihood(jlmm) atol=0.001
        end

        @testset "Gaussian" begin
            # TODO: remove upon next MixedModels.jl release
            # TODO: dispersion check
            # dat = dataset(:verbagg);
            # jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)+(1|item)),
            #             dat, Bernoulli());
            #
            # jm = (jlmm, dat);
            # fit!(jlmm, fast=true)
            # @rput jm;
            # # @test_warn Regex(".*categorical.*") @rput jm;
            # # note the really high tolerances
            # @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.1
            # @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.5
        end

        @testset "InverseGaussian" begin
        end

        @testset "Gamma" begin
        end

        @testset "contrasts" begin
            jlmm = fit(MixedModel, @formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)+(1|item)),
                        dat, Bernoulli(), contrasts=Dict(:gender => EffectsCoding(),
                                                         :btypes => EffectsCoding()));
            jm = Tuple([jlmm, dat]);
            @rput jm;
            @test fixef(jlmm) ≈ rcopy(R"fixef(jm)");
        end

    end
end
