using RCall, MixedModels, Test
using StatsBase: zscore
using StatsModels: SeqDiffCoding
using Tables: columntable

import JellyMe4: _set_lmer, _set_afex_installed

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
    }
    if(!require(afex)){
        .libPaths("/tmp")
        lib <- .libPaths()[1L]
        install.packages("afex",repos="https://cloud.r-project.org", libs=lib)
    }
    """)

    # this is available in MixedModels.dataset(:sleepstudy) but with different
    # capitalization than in R
    sleepstudy = rcopy(R"sleepstudy")

    @testset "get lmerMod" begin
        ### from R ###

        jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + (1 + Days|Subject)),sleepstudy), REML=false)
        rlmm = rcopy(R"m <- lme4::lmer(Reaction ~ 1 + Days + (1 + Days|Subject),sleepstudy,REML=FALSE)")

        @test jlmm.θ ≈ rlmm.θ atol=0.001
        @test objective(jlmm) ≈ objective(rlmm) atol=0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001

        jlmm = refit!(jlmm, REML=true)
        rlmm = rcopy(R"update(m, REML=TRUE)")

        @test jlmm.θ ≈ rlmm.θ atol=0.001
        @test objective(jlmm) ≈ objective(rlmm) atol=0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001

        @testset "merModLmerTest" begin
            jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + (1 + Days|Subject)),sleepstudy), REML=false)
            rlmm = rcopy(R"m <- lmerTest::lmer(Reaction ~ 1 + Days + (1 + Days|Subject),sleepstudy,REML=FALSE)")

            @test jlmm.θ ≈ rlmm.θ atol=0.001
            @test objective(jlmm) ≈ objective(rlmm) atol=0.001
            @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001

            jlmm = refit!(jlmm, REML=true)
            rlmm = rcopy(R"update(m, REML=TRUE)")

            @test jlmm.θ ≈ rlmm.θ atol=0.001
            @test objective(jlmm) ≈ objective(rlmm) atol=0.001
            @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001
        end

        # single scalar RE is different because RCall converts the single-entry θ
        # to a scalar value, which we have to correct
        # of course, this isn't a problem in the other direction, because
        # the scalar-vector distinction for θ is missing in R
        @testset "scalar RE" begin
            jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + (1|Subject)),sleepstudy), REML=false)
            rlmm = rcopy(R"m <- lme4::lmer(Reaction ~ 1 + Days + (1|Subject),sleepstudy,REML=FALSE)")

            @test jlmm.θ ≈ rlmm.θ atol=0.001
            @test objective(jlmm) ≈ objective(rlmm) atol=0.001
            @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001
        end

        @testset "contrasts" begin
            reval("""
            data(cake)
            cake\$rr <- with(cake, replicate:recipe)
            """);
            rlmm = rcopy(R"fm1 <- lme4::lmer(angle ~ recipe * temperature + (1|rr), cake, REML= FALSE)");
            @test fixef(rlmm) ≈  rcopy(R"fixef(fm1)");
            # rlmm = rcopy(R"""fm1 <- lme4::lmer(angle ~ recipe * temperature + (1|rr), cake, REML= FALSE,
            #                              contrasts=list(temperature=contr.helmert))""");
            # @test fixef(rlmm) ≈  rcopy(R"fixef(fm1)");
        end

        @testset "double-bar" begin
            # double bar works because lme4 treats this internally as a rewrite rule into distinct terms
            # the distinct terms are handled the same way in both languages, so everything is fine.
            # printing differs a lot between languages, but that's life


            reval("""
            machines <- as.data.frame(nlme::Machines)
            mach <- lme4::lmer(score ~ Machine + (Machine || Worker), machines, REML=FALSE)
            """)
            machines = rcopy(R"machines")
            @test_throws ArgumentError rlmm = rcopy(R"mach")

            # note that is the lme4 definition of double || -- it really is just a convenience wrapper for
            # splitting terms up this way!
            # jlmm = fit(MixedModel, @formula(score ~ 1 +  Machine + (1|Worker) + (fulldummy(Machine)|Worker)), machines)
            # as a cheat for comparing the covariance matrices, we use PCA
            # @test only(rlmm.rePCA) ≈ only(jlmm.rePCA) atol=0.05

            rlmm = rcopy(R"m <- lme4::lmer(Reaction ~ 1 + Days + (1 + Days||Subject),sleepstudy,REML=FALSE)")
            jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + zerocorr(1 + Days|Subject)),sleepstudy), REML=false)
            # as a cheat for comparing the covariance matrices, we use PCA
            @test only(rlmm.rePCA) ≈ only(jlmm.rePCA) atol=0.05
        end

        @testset "dummy" begin
            # TODO
        end
    end

    @testset "put lmerMod" begin
        ### from Julia ###
        jlmm = LMM(@formula(Reaction ~ 1 + Days + (1 + Days|Subject)),sleepstudy)
        jm = (jlmm, sleepstudy)
        # unfitted model
        @test_throws ArgumentError @rput jm
        fit!(jlmm, REML=true)
        @rput jm
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
        @test rcopy(R"REMLcrit(jm)") ≈ objective(jlmm)

        refit!(jlmm, REML=false)
        jm = (jlmm, sleepstudy)
        @rput jm
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
        @test rcopy(R"deviance(jm)") ≈ objective(jlmm)

        @testset "columntable" begin
            jm = (jlmm, columntable(sleepstudy));
            @rput jm;
        end

        @testset "transformations" begin
            sleepstudy[!,:Days2] = sleepstudy.Days .+ 1
            @rput sleepstudy;
            R"m <- lme4::lmer(log10(Reaction) ~ 1 + log(Days2) + (1 + log(Days2)|Subject),sleepstudy,REML=FALSE)"
            jlmm = fit!(LMM(@formula(log10(Reaction) ~ 1 + log(Days2) + (1 + log(Days2)|Subject)),sleepstudy), REML=false)
            jm = (jlmm, sleepstudy)
            @rput jm;
            @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
            @test rcopy(R"deviance(jm)") ≈ objective(jlmm)

            jlmm = fit!(LMM(@formula(Reaction ~ 1 + round(Days) + (1|Subject)),sleepstudy), REML=false)
            jm = (jlmm, sleepstudy);

            @test_throws ArgumentError (@rput jm)
        end

        @testset "contrasts" begin
            cake = rcopy(reval("""
            data(cake)
            cake\$rr <- with(cake, replicate:recipe)
            cake
            """));
            jlmm = fit(MixedModel, @formula(angle ~ recipe * temperature + (1|rr)),
                       cake, REML=false, contrasts=Dict(:temperature => SeqDiffCoding()));
            jm = (jlmm, cake);
            @rput jm;
            @test fixef(jlmm) ≈  rcopy(R"fixef(jm)");
        end


        @testset "fulldummy" begin
            machines = rcopy(R"as.data.frame(nlme::Machines)")
            jlmm = fit(MixedModel, @formula(score ~ 1 +  Machine + (1 + fulldummy(Machine)|Worker)), machines)
            rlmm = (jlmm, machines)
            @rput rlmm
            rlmmrepca = rcopy(R"summary(rePCA(rlmm))$Worker$importance[3,]")
            @test  rlmmrepca ≈ only(MixedModels.rePCA(jlmm; corr=false)) atol=0.05
        end

        @testset "zerocorr" begin
            # TODO: test fulldummy within a zerocorr when zerocorr is better supported


            @testset "lme4" begin
                _set_lmer("lme4::lmer")
                _set_afex_installed(false)

                jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + zerocorr(1 + Days|Subject)),sleepstudy), REML=false)
                rlmm = (jlmm, sleepstudy)
                @rput rlmm
                @test rcopy(R"""!is(rlmm,"merModLmerTest")""")
                @test only(ranef(jlmm))' ≈ Matrix(rcopy(R"ranef(rlmm)$Subject"))
                @test fixef(jlmm) ≈ rcopy(R"fixef(rlmm)")
                @test vcov(jlmm) ≈ rcopy(R"as.matrix(vcov(rlmm))")
            end


            @testset "afex" begin

                machines = rcopy(R"as.data.frame(nlme::Machines)")

                jlmm = fit(MixedModel, @formula(score ~ 1 +  Machine + zerocorr(0+Machine|Worker)), machines)

                @testset "afex pre-enabled" begin
                    _set_lmer("afex::lmer_alt")
                    rlmm = (jlmm, machines)
                    @rput rlmm
                    @test only(ranef(jlmm))' ≈ Matrix(rcopy(R"ranef(rlmm)$Worker"))
                    @test fixef(jlmm) ≈ rcopy(R"fixef(rlmm)")
                    @test vcov(jlmm) ≈ rcopy(R"as.matrix(vcov(rlmm))")
                end

                @testset "afex auto-enabled" begin
                    _set_lmer("lme4::lmer")
                    _set_afex_installed(true)
                    rlmm = (jlmm, machines)
                    # match_mode needs to specified in case there are further
                    # warnings from R/RCall
                    @test_logs (:info, r"afex::lmer_alt") match_mode=:any @rput rlmm
                    @test only(ranef(jlmm))' ≈ Matrix(rcopy(R"ranef(rlmm)$Worker"))
                    @test fixef(jlmm) ≈ rcopy(R"fixef(rlmm)")
                    @test vcov(jlmm) ≈ rcopy(R"as.matrix(vcov(rlmm))")
                end

                @testset "afex unavailable" begin
                    _set_lmer("lme4::lmer")
                    _set_afex_installed(false)
                    rlmm = (jlmm, machines)
                    @test_throws ArgumentError @rput rlmm
                end
            end

        end
    end
end
