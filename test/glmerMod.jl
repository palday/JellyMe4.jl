using RCall, MixedModels, Test
using StatsBase: zscore
using Tables: columntable
const LMM = LinearMixedModel
const GLMM = GeneralizedLinearMixedModel

logistic(x)  = 1 / (1 + exp(-x))

@testset "glmerMod" begin
    reval("""
    if(!require(lme4)){
        # use tmp for tests if lme4 isn't available
        .libPaths("/tmp")
        lib <- .libPaths()[1L]
        install.packages("lme4",repos="https://cloud.r-project.org", libs=lib)
        library(lme4)
    }""")

    ### from R ###
    @testset "get glmerMod" begin
        # from lme4 ?cbpp
        R"""
        ## response as a matrix
        m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                      family = binomial, data = cbpp, nAGQ=0)
         ## response as a vector of probabilities and usage of argument "weights"
         m1p <- glmer(incidence / size ~ period + (1 | herd), weights = size,
                      family = binomial, data = cbpp)
        """

    end

    ### from Julia ###
    @testset "put glmerMod" begin
        @testset "Bernoulli" begin
            # NOTE: For Bernoulli, some of the tests are doubled in order to
            #       check nAQG and fast options
            # TODO: remove upon next MixedModels.jl release
            dat = dataset(:verbagg);

            @testset "unfitted model" begin
                # we max out at 1 scalar RE for the nAGQ tests
                jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)),
                            dat, Bernoulli());
                jm = Tuple([jlmm, dat]);
                # unfitted model
                @test_throws ArgumentError @rput jm
            end

            @testset "nAGQ" begin
                # we max out at 1 scalar RE for the nAGQ tests
                jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)),
                            dat, Bernoulli());
                fit!(jlmm; fast=false, nAGQ=9)
                jm = (jlmm, dat);
                @rput jm;
                # @test_warn Regex(".*categorical.*") @rput jm;

                @test  rcopy(R"""jm@devcomp$dims["nAGQ"]""") == jlmm.optsum.nAGQ
                # note the really high tolerances
                @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.1
                @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.1
            end

            @testset "Laplace" begin
                # we max out at 1 scalar RE for the nAGQ tests
                jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)),
                            dat, Bernoulli());
                fit!(jlmm);
                jm = Tuple([jlmm, dat])
                @rput jm;

                # note the really high tolerances
                @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.1
                @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.5;
            end

            @testset "columntable" begin
                jlmm = GLMM(@formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)),
                            dat, Bernoulli());
                fit!(jlmm);
                jm = Tuple([jlmm, columntable(dat)]);
                @rput jm;
            end

        end

        @testset "Binomial" begin
            # TODO: remove upon next MixedModels.jl release
            dat  = dataset(:cbpp);
            dat.rate = dat.incid ./ dat.hsz;
            jlmm = fit(MixedModel, @formula(rate ~ 1 + period + (1|herd)),
                      dat, Binomial(); wts=float(dat.hsz), fast=true);

            jm = (jlmm, dat);
            @rput jm;
            # @test_warn Regex(".*categorical.*") @rput jm;
            @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol=0.001
            @test_broken rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol=0.001
            # this is where the weights have their biggest impact
            @test stderror(jlmm) ≈ rcopy(R"""coef(summary(jm))[,"Std. Error"]""") atol=0.001
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
            dat = dataset(:verbagg);

            jlmm = fit(MixedModel, @formula(r2 ~ 1+anger+gender+btype+situ+(1|subj)+(1|item)),
                        dat, Bernoulli(), contrasts=Dict(:gender => EffectsCoding(),
                                                         :btypes => EffectsCoding()));
            jm = Tuple([jlmm, dat]);
            @rput jm;
            @test fixef(jlmm) ≈ rcopy(R"fixef(jm)");
        end
    end
end
