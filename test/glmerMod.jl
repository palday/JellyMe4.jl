using RCall, MixedModels, Test
using StatsBase: zscore
using Tables: columntable
using GLM
import GLM: Link
const LMM = LinearMixedModel
const GLMM = GeneralizedLinearMixedModel

logistic(x) = 1 / (1 + exp(-x))

@testset "glmerMod" begin
    reval("""
    if(!("lme4" %in% installed.packages())){
        # use tmp for tests if lme4 isn't available
        .libPaths("/tmp")
        lib <- .libPaths()[1L]
        install.packages("lme4",repos="https://cloud.r-project.org", libs=lib)
    }
    if(!("afex" %in% installed.packages())){
        .libPaths("/tmp")
        lib <- .libPaths()[1L]
        install.packages("afex",repos="https://cloud.r-project.org", libs=lib)
    }
    """)

    ### from R ###
    @testset "get glmerMod" begin
        @testset "Binomial" begin
            reval(raw"""
            attach(cbpp)
            cbpp$rate <- with(cbpp, incidence/size)
            rlmm <- glmer(rate ~ period + (1 | herd), weights = size,
                          family = binomial, data = cbpp)
            """)
            rlmm = rcopy(R"rlmm")
            @test rcopy(R"fitted(rlmm)") ≈ fitted(rlmm) atol = 0.001
            @test_broken rcopy(R"-2 * logLik(rlmm)") ≈ deviance(rlmm) atol = 0.001
            # this is where the weights have their biggest impact
            @test stderror(rlmm) ≈ rcopy(R"""coef(summary(rlmm))[,"Std. Error"]""") atol = 0.1
            @test rcopy(R"logLik(rlmm)") ≈ loglikelihood(rlmm) atol = 0.001

            @testset "alternative response specifications" begin
                # R"""
                # # from lme4 ?cbpp
                # ## response as a matrix
                # m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                #               family = binomial, data = cbpp, nAGQ=0)
                #  ## response as a vector of probabilities and usage of argument "weights"
                #  m1p <- glmer(incidence / size ~ period + (1 | herd), weights = size,
                #               family = binomial, data = cbpp)
                # """

                @testset "cbind" begin end

                @testset "proportion computed in line" begin end
            end

            @testset "contrasts" begin
                reval(raw"""
                contrasts(cbpp$period) <- contr.helmert(levels(cbpp$period))
                """)
                rlmm = rcopy(R"fm1 <- glmer(rate ~ period + (1 | herd), weights = size,family = binomial, data = cbpp)")
                @test fixef(rlmm) ≈ rcopy(R"fixef(fm1)")
            end

            @testset "Bernoulli" begin
                @testset "nAGQ and scalar RE" begin
                    reval("""
                    rlmm <- glmer(r2 ~ Anger + Gender + btype + situ + (1|id),
                                        family = binomial, data=VerbAgg, nAGQ=4)
                    """)

                    rlmm = rcopy(R"rlmm")
                    @test rcopy(R"fitted(rlmm)") ≈ fitted(rlmm) atol = 0.001
                    @test rcopy(R"-2 * logLik(rlmm)") ≈ deviance(rlmm) atol = 0.001
                    # this is where the weights have their biggest impact
                    @test stderror(rlmm) ≈ rcopy(R"""coef(summary(rlmm))[,"Std. Error"]""") atol = 0.05
                    @test_broken rcopy(R"logLik(rlmm)") ≈ loglikelihood(rlmm) atol = 0.001
                end

                @testset "Laplace and probit link" begin
                    # TODO add tests for each link
                    reval("""
                    rlmm <- glmer(r2 ~ Anger + Gender + btype + situ + (1|id),
                                        family = binomial(link="probit"), data=VerbAgg)
                    """)

                    rlmm = rcopy(R"rlmm")
                    @test Link(rlmm.resp) isa GLM.ProbitLink
                    @test rcopy(R"fitted(rlmm)") ≈ fitted(rlmm) atol = 0.001
                    @test_broken rcopy(R"-2 * logLik(rlmm)") ≈ deviance(rlmm) atol = 0.001
                    # this is where the weights have their biggest impact
                    @test stderror(rlmm) ≈ rcopy(R"""coef(summary(rlmm))[,"Std. Error"]""") atol = 0.05
                    @test rcopy(R"logLik(rlmm)") ≈ loglikelihood(rlmm) atol = 0.001
                end
            end
        end

        @testset "Poisson and fast fit" begin
            reval("""
            # from the lme4 docs
            form <- TICKS ~ YEAR + HEIGHT + (1|BROOD) + (1|INDEX) + (1|LOCATION)
            rlmm  <- glmer(form, family="poisson",data=grouseticks, nAGQ=0)
            """)
            rlmm = rcopy(R"rlmm")
            @test rcopy(R"fitted(rlmm)") ≈ fitted(rlmm) atol = 0.001
            @test_broken rcopy(R"-2 * logLik(rlmm)") ≈ deviance(rlmm) atol = 0.001
            # 7.5% difference isn't great....
            @test stderror(rlmm) ≈ rcopy(R"""coef(summary(rlmm))[,"Std. Error"]""") rtol = 0.075
            @test rcopy(R"logLik(rlmm)") ≈ loglikelihood(rlmm) atol = 0.05
        end
    end

    ### from Julia ###
    @testset "put glmerMod" begin
        @testset "Bernoulli" begin
            # NOTE: For Bernoulli, some of the tests are doubled in order to
            #       check nAQG and fast options
            # TODO: remove upon next MixedModels.jl release
            dat = dataset(:verbagg)

            @testset "unfitted model" begin
                # we max out at 1 scalar RE for the nAGQ tests
                jlmm = GLMM(@formula(r2 ~ 1 + anger + gender + btype + situ + (1 | subj)),
                            dat, Bernoulli())
                jm = Tuple([jlmm, dat])
                # unfitted model
                @test_throws ArgumentError @rput jm
            end

            @testset "nAGQ" begin
                # we max out at 1 scalar RE for the nAGQ tests
                jlmm = GLMM(@formula(r2 ~ 1 + anger + gender + btype + situ + (1 | subj)),
                            dat, Bernoulli())
                fit!(jlmm; fast=false, nAGQ=9, progress=false)
                jm = (jlmm, dat)
                @rput jm
                # @test_warn Regex(".*categorical.*") @rput jm;

                @test rcopy(R"""jm@devcomp$dims["nAGQ"]""") == jlmm.optsum.nAGQ
                # note the really high tolerances
                @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol = 0.1
                @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol = 0.1
            end

            @testset "Laplace" begin
                # we max out at 1 scalar RE for the nAGQ tests
                jlmm = GLMM(@formula(r2 ~ 1 + anger + gender + btype + situ + (1 | subj)),
                            dat, Bernoulli())
                fit!(jlmm; progress=false)
                jm = Tuple([jlmm, dat])
                @rput jm

                # note the really high tolerances
                @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol = 0.1
                @test rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol = 0.5
            end

            @testset "columntable" begin
                jlmm = GLMM(@formula(r2 ~ 1 + anger + gender + btype + situ + (1 | subj)),
                            dat, Bernoulli())
                fit!(jlmm; progress=false)
                jm = Tuple([jlmm, columntable(dat)])
                @rput jm
            end
        end

        @testset "Binomial" begin
            # TODO: remove upon next MixedModels.jl release
            dat = dataset(:cbpp)
            dat.rate = dat.incid ./ dat.hsz
            jlmm = fit(MixedModel, @formula(rate ~ 1 + period + (1 | herd)),
                       dat, Binomial(); wts=float(dat.hsz), fast=true, progress=false)

            jm = (jlmm, dat)
            @rput jm
            # @test_warn Regex(".*categorical.*") @rput jm;
            @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol = 0.001
            @test_broken rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol = 0.001
            # this is where the weights have their biggest impact
            @test stderror(jlmm) ≈ rcopy(R"""coef(summary(jm))[,"Std. Error"]""") atol = 0.1
            @test rcopy(R"logLik(jm)") ≈ loglikelihood(jlmm) atol = 0.001
        end

        @testset "Poisson and fast fit" begin
            # TODO: remove upon next MixedModels.jl release
            center(v::AbstractVector) = v .- (sum(v) / length(v))
            dat = dataset(:grouseticks)
            dat.ch = center(dat.height)
            jlmm = GLMM(@formula(ticks ~ 1 + year + ch + (1 | index) + (1 | brood) +
                                         (1 | location)),
                        dat, Poisson())

            jm = (jlmm, dat)
            fit!(jlmm; fast=true, progress=false) # problems with this one in fast=false
            @rput jm
            # @test_warn Regex(".*categorical.*") @rput jm;
            @test rcopy(R"""jm@devcomp$dims["nAGQ"]""") == 0
            # note the really high tolerances
            @test rcopy(R"fitted(jm)") ≈ fitted(jlmm) atol = 0.001
            @test_broken rcopy(R"-2 * logLik(jm)") ≈ deviance(jlmm) atol = 0.5
            @test rcopy(R"logLik(jm)") ≈ loglikelihood(jlmm) atol = 0.001
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

        @testset "InverseGaussian" begin end

        @testset "Gamma" begin end

        @testset "contrasts" begin
            dat = dataset(:verbagg)

            jlmm = fit(MixedModel,
                       @formula(r2 ~ 1 + anger + gender + btype + situ + (1 | subj) +
                                     (1 | item)),
                       dat, Bernoulli(); progress=false,
                       contrasts=Dict(:gender => EffectsCoding(),
                                      :btypes => EffectsCoding()))
            jm = (jlmm, dat)
            @rput jm
            @test fixef(jlmm) ≈ rcopy(R"fixef(jm)")

          @testset "asinh transformation" begin
            dat = dataset(:verbagg)

            jlmm = fit(MixedModel,
                       @formula(r2 ~ 1 + asinh(anger) + gender + btype + situ + (1 | subj) +
                                     (1 | item)),
                       dat, Bernoulli())
            jm = (jlmm, dat)
            @rput jm
            @test fixef(jlmm) ≈ rcopy(R"fixef(jm)")
        end
    end
end
