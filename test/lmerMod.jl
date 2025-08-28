# the associated Julia dataset has different capitalization
# capitalization than in R
sleepstudy = rcopy(R"lme4::sleepstudy")

kb07 = dataset(:kb07)

@testset "get lmerMod" begin
    ### from R ###

    jlmm = lmm(@formula(Reaction ~ 1 + Days + (1 + Days | Subject)), sleepstudy;
               REML=false, progress=false)
    rlmm = rcopy(R"m <- lme4::lmer(Reaction ~ 1 + Days + (1 + Days|Subject),lme4::sleepstudy,REML=FALSE)")

    @test jlmm.θ ≈ rlmm.θ atol = 0.001
    @test objective(jlmm) ≈ objective(rlmm) atol = 0.001
    @test fixef(jlmm) ≈ fixef(rlmm) atol = 0.001

    jlmm = refit!(jlmm; REML=true, progress=false)
    rlmm = rcopy(R"update(m, REML=TRUE)")

    @test jlmm.θ ≈ rlmm.θ atol = 0.001
    @test objective(jlmm) ≈ objective(rlmm) atol = 0.001
    @test fixef(jlmm) ≈ fixef(rlmm) atol = 0.001

    @testset "merModLmerTest" begin
        jlmm = lmm(@formula(Reaction ~ 1 + Days + (1 + Days | Subject)),
                   sleepstudy; REML=false, progress=false)
        rlmm = rcopy(R"m <- lmerTest::lmer(Reaction ~ 1 + Days + (1 + Days|Subject),sleepstudy,REML=FALSE)")

        @test jlmm.θ ≈ rlmm.θ atol = 0.001
        @test objective(jlmm) ≈ objective(rlmm) atol = 0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol = 0.001

        jlmm = refit!(jlmm; REML=true, progress=false)
        rlmm = rcopy(R"update(m, REML=TRUE)")

        @test jlmm.θ ≈ rlmm.θ atol = 0.001
        @test objective(jlmm) ≈ objective(rlmm) atol = 0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol = 0.001
    end

    # single scalar RE is different because RCall converts the single-entry θ
    # to a scalar value, which we have to correct
    # of course, this isn't a problem in the other direction, because
    # the scalar-vector distinction for θ is missing in R
    @testset "scalar RE" begin
        jlmm = lmm(@formula(Reaction ~ 1 + Days + (1 | Subject)), sleepstudy;
                   REML=false, progress=false)
        rlmm = rcopy(R"m <- lme4::lmer(Reaction ~ 1 + Days + (1|Subject),sleepstudy,REML=FALSE)")

        @test jlmm.θ ≈ rlmm.θ atol = 0.001
        @test objective(jlmm) ≈ objective(rlmm) atol = 0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol = 0.001
    end

    @testset "sorting by n BLUPs vs. n groups" begin
        @rput kb07
        # ignore the boundary fit
        rlmm = @suppress rcopy(R"rlmm <- lmer(rt_trunc ~ 1 + (1|subj)+(1+load+prec+spkr|item), $(kb07), REML=FALSE)")
        @test rcopy(R"fitted(rlmm)") ≈ fitted(rlmm)
        @test rcopy(R"deviance(rlmm)") ≈ objective(rlmm)
        @test all(isapprox.(collect(VarCorr(rlmm).σρ.item.σ),
                            rcopy(R"""attr(VarCorr(rlmm)[["item"]], "stddev")""")))
        @test all(isapprox.(collect(VarCorr(rlmm).σρ.subj.σ),
                            rcopy(R"""attr(VarCorr(rlmm)[["subj"]], "stddev")""")))
    end

    @testset "contrasts" begin
        reval("""
        cake <- lme4::cake
        cake\$rr <- with(cake, replicate:recipe)
        """)
        rlmm = @suppress rcopy(R"fm1 <- lme4::lmer(angle ~ recipe * temperature + (1|rr), cake, REML= FALSE)")
        @test fixef(rlmm) ≈ rcopy(R"fixef(fm1)")
        # rlmm = rcopy(R"""fm1 <- lme4::lmer(angle ~ recipe * temperature + (1|rr), cake, REML= FALSE,
        #                              contrasts=list(temperature=contr.helmert))""");
        # @test fixef(rlmm) ≈  rcopy(R"fixef(fm1)");
    end

    @testset "caret" begin
        reval("""
        cake <- lme4::cake
        cake\$rr <- with(cake, replicate:recipe)
        """)
        rlmm = @suppress rcopy(R"fm1 <- lme4::lmer(angle ~ (recipe + temperature)^2 + (1|rr), cake, REML= FALSE)")
        @test fixef(rlmm) ≈ rcopy(R"fixef(fm1)")
    end

    @testset "double-bar" begin
        # double bar works because lme4 treats this internally as a rewrite rule into distinct terms
        # the distinct terms are handled the same way in both languages, so everything is fine.
        # printing differs a lot between languages, but that's life

        @suppress reval("""
        machines <- as.data.frame(nlme::Machines)
        mach <- lme4::lmer(score ~ Machine + (Machine || Worker), machines, REML=FALSE)
        """)
        machines = rcopy(R"machines")
        @suppress @test_throws ArgumentError rcopy(R"mach")

        # note that is the lme4 definition of double || -- it really is just a convenience wrapper for
        # splitting terms up this way!
        # jlmm = fit(MixedModel, @formula(score ~ 1 +  Machine + (1|Worker) + (fulldummy(Machine)|Worker)), machines)
        # as a cheat for comparing the covariance matrices, we use PCA
        # @test only(rlmm.rePCA) ≈ only(jlmm.rePCA) atol=0.05

        rlmm = rcopy(R"m <- lme4::lmer(Reaction ~ 1 + Days + (1 + Days||Subject),sleepstudy,REML=FALSE)")
        jlmm = lmm(@formula(Reaction ~ 1 + Days + zerocorr(1 + Days | Subject)),
                   sleepstudy; REML=false, progress=false)
        # as a cheat for comparing the covariance matrices, we use PCA
        @test only(rlmm.rePCA) ≈ only(jlmm.rePCA) atol = 0.05
    end

    @testset "nested grouping" begin
        pastes = DataFrame(dataset(:pastes))
        @rput pastes
        jlmm = lmm(@formula(strength ~ 1 + (1 | batch / cask)), pastes;
                   REML=false, progress=false)
        rlmm = @suppress rcopy(R"lme4::lmer(strength ~ 1 + (1 | batch / cask), pastes,REML=FALSE)")

        @test jlmm.θ ≈ rlmm.θ atol = 0.001
        @test objective(jlmm) ≈ objective(rlmm) atol = 0.001
        @test fixef(jlmm) ≈ fixef(rlmm) atol = 0.001
    end

    @testset "dummy" begin
        # TODO
    end
end

@testset "put lmerMod" begin
    ### from Julia ###
    jlmm = LinearMixedModel(@formula(Reaction ~ 1 + Days + (1 + Days | Subject)),
                            sleepstudy)
    jm = (jlmm, sleepstudy)
    # unfitted model
    @test_throws ArgumentError @rput(jm)
    fit!(jlmm; REML=true, progress=false)
    @rput jm
    @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
    @test rcopy(R"REMLcrit(jm)") ≈ objective(jlmm)

    refit!(jlmm; REML=false, progress=false)
    jm = (jlmm, sleepstudy)
    @rput jm
    @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
    @test rcopy(R"deviance(jm)") ≈ objective(jlmm)

    @testset "columntable" begin
        jm = (jlmm, columntable(sleepstudy))
        @rput jm
        @test true
    end

    @testset "transformations" begin
        sleepstudy[!, :Days2] = sleepstudy.Days .+ 1
        @rput sleepstudy
        @suppress R"m <- lme4::lmer(log10(Reaction) ~ 1 + log(Days2) + (1 + log(Days2)|Subject),sleepstudy,REML=FALSE)"
        jlmm = lmm(@formula(log10(Reaction) ~ 1 + log(Days2) +
                                              (1 + log(Days2) | Subject)),
                   sleepstudy; REML=false, progress=false)
        jm = (jlmm, sleepstudy)
        @rput jm
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
        @test rcopy(R"deviance(jm)") ≈ objective(jlmm)

        jlmm = lmm(@formula(Reaction ~ 1 + round(Days) + (1 | Subject)),
                   sleepstudy; REML=false, progress=false)
        jm = (jlmm, sleepstudy)

        @test_throws ArgumentError @rput(jm)
    end

    @testset "sorting by n BLUPs vs. n groups" begin
        jlmm = lmm(@formula(rt_trunc ~ 1 + (1 | subj) + (1 + load + prec + spkr | item)),
                   kb07; progress=false)
        jm = (jlmm, kb07)
        # suppress singular warning
        @suppress @rput jm
        @test rcopy(R"fitted(jm)") ≈ fitted(jlmm)
        @test rcopy(R"deviance(jm)") ≈ objective(jlmm)
        @test all(isapprox.(collect(VarCorr(jlmm).σρ.item.σ),
                            rcopy(R"""attr(VarCorr(jm)[["item"]], "stddev")""")))
        @test all(isapprox.(collect(VarCorr(jlmm).σρ.subj.σ),
                            rcopy(R"""attr(VarCorr(jm)[["subj"]], "stddev")""")))
    end

    @testset "contrasts" begin
        cake = rcopy(reval("""
        data(cake)
        cake\$rr <- with(cake, replicate:recipe)
        cake
        """))
        jlmm = lmm(@formula(angle ~ recipe * temperature + (1 | rr)),
                   cake; REML=false, progress=false,
                   contrasts=Dict(:temperature => SeqDiffCoding()))
        jm = (jlmm, cake)
        # suppress contrast info
        @suppress @rput jm
        @test fixef(jlmm) ≈ rcopy(R"fixef(jm)")
    end

    @testset "fulldummy" begin
        machines = rcopy(R"as.data.frame(nlme::Machines)")
        jlmm = lmm(@formula(score ~ 1 + Machine + (1 + fulldummy(Machine) | Worker)),
                   machines; progress=false)
        rlmm = (jlmm, machines)
        @rput rlmm
        rlmmrepca = rcopy(R"summary(rePCA(rlmm))$Worker$importance[3,]")
        @test rlmmrepca ≈ only(MixedModels.rePCA(jlmm; corr=false)) atol = 0.05
    end

    @testset "zerocorr" begin
        # TODO: test fulldummy within a zerocorr when zerocorr is better supported

        @testset "lme4" begin
            _set_lmer("lme4::lmer")
            _set_afex_installed(false)

            jlmm = lmm(@formula(Reaction ~ 1 + Days + zerocorr(1 + Days | Subject)),
                       sleepstudy; REML=false, progress=false)
            rlmm = (jlmm, sleepstudy)
            @rput rlmm
            @test rcopy(R"""!is(rlmm,"merModLmerTest")""")
            @test only(ranef(jlmm))' ≈ Matrix(rcopy(R"ranef(rlmm)$Subject"))
            @test fixef(jlmm) ≈ rcopy(R"fixef(rlmm)")
            @test vcov(jlmm) ≈ rcopy(R"as.matrix(vcov(rlmm))")
        end

        @testset "afex" begin
            machines = rcopy(R"as.data.frame(nlme::Machines)")

            jlmm = lmm(@formula(score ~ 1 + Machine + zerocorr(0 + Machine | Worker)),
                       machines; progress=false)

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
                @test_logs (:info, r"afex::lmer_alt") match_mode = :any @rput(rlmm)
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

    @testset "nested grouping" begin
        pastes = DataFrame(dataset(:pastes))
        @rput pastes
        jlmm = lmm(@formula(strength ~ 1 + (1 | batch / cask)), pastes;
                   REML=false, progress=false)
        rlmm = (jlmm, pastes)
        @rput rlmm

        @test jlmm.θ ≈ rcopy(R"""getME(rlmm, "theta")""") atol = 0.001
        @test objective(jlmm) ≈ rcopy(R"deviance(rlmm)") atol = 0.001
        @test only(fixef(jlmm)) ≈ rcopy(R"fixef(rlmm)") atol = 0.001
    end
end

@testset "lmerControl warnings" begin
    slp = DataFrame(MixedModels.dataset("sleepstudy"))
    slp[!, :obs] .= 1:nrow(slp)
    fm1 = lmm(@formula(reaction ~ 1 + days + (1 | obs)),
              slp; progress=false,
              contrasts=Dict(:obs => Grouping()))
    rfm1 = (fm1, slp)
    warning = r"""Warning: number of levels of each grouping factor must be < number of observations \(problems: obs\)
                  Warning: number of observations \(=180\) <= number of random effects \(=180\) for term \(1 \| obs\); the random-effects parameters and the residual variance \(or scale parameter\) are probably unidentifiable
                  Warning: number of observations \(=180\) <= rank\(Z\) \(=180\); the random-effects parameters and the residual variance \(or scale parameter\) are probably unidentifiable
                  """
    @test_logs (:warn, warning) @rput rfm1
end
