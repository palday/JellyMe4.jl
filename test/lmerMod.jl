using RCall, MixedModels, Test
using StatsBase: zscore
import StatsModels: SeqDiffCoding
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

        @testset "contrasts" begin
            reval("""
            data(cake)
            cake\$rr <- with(cake, replicate:recipe)
            """);
            rlmm = rcopy(R"fm1 <- lmer(angle ~ recipe * temperature + (1|rr), cake, REML= FALSE)");
            @test fixef(rlmm) ≈  rcopy(R"fixef(fm1)");
            # rlmm = rcopy(R"""fm1 <- lmer(angle ~ recipe * temperature + (1|rr), cake, REML= FALSE,
            #                              contrasts=list(temperature=contr.helmert))""");
            # @test fixef(rlmm) ≈  rcopy(R"fixef(fm1)");
        end
        
        @testset "double-bar" begin
            # double bar works because lme4 treats this internally as a rewrite rule into distinct terms
            # the distinct terms are handled the same way in both languages, so everything is fine.
            # printing differs a lot between languages, but that's life
            
            
            reval("""
            machines <- as.data.frame(nlme::Machines)
            mach <- lmer(score ~ Machine + (Machine || Worker), machines, REML=FALSE)
            """)
            machines = rcopy(R"machines")
            rlmm = rcopy(R"mach")
            jlmm = fit(MixedModel, @formula(score ~ 1 +  Machine + (1|Worker) + (0+Machine|Worker)), machines)
            # as a cheat for comparing the covariance matrices, we use packages
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

        fit!(jlmm, REML=false)
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
            R"m <- lmer(log10(Reaction) ~ 1 + log(Days2) + (1 + log(Days2)|Subject),sleepstudy,REML=FALSE)"
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
            rlmm = R"lmer(score ~ 1 +  Machine + (1 + dummy(Machine, levels(Machine))|Worker), machines)"
            rlmmrepca = rcopy(R"summary(rePCA($rlmm))$Worker$importance[3,]")
            @test  rlmmrepca ≈ only(jlmm.rePCA) atol=0.05
        end
                
        @testset "zerocorr" begin    
            jlmm = fit!(LMM(@formula(Reaction ~ 1 + Days + zerocorr(1 + Days|Subject)),sleepstudy), REML=false)
            rlmm = rcopy(R"m <- lmer(Reaction ~ 1 + Days + (1 + Days||Subject),sleepstudy,REML=FALSE)")

            @test jlmm.θ ≈ rlmm.θ atol=0.001
            @test objective(jlmm) ≈ objective(rlmm) atol=0.001
            @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001

            jlmm = fit!(jlmm, REML=true)
            rlmm = rcopy(R"update(m, REML=TRUE)")

            @test jlmm.θ ≈ rlmm.θ atol=0.001
            @test objective(jlmm) ≈ objective(rlmm) atol=0.001
            @test fixef(jlmm) ≈ fixef(rlmm) atol=0.001
            # TODO: test fulldummy within a zerocorr when zerocorr is better supported

            # machines = rcopy(R"as.data.frame(nlme::Machines)")
            # jlmm = fit(MixedModel, @formula(score ~ 1 +  Machine + zerocorr(1+Machine|Worker)), machines)
            
            
            # rlmm = rcopy(R"mach")
            # jlmm = fit(MixedModel, @formula(score ~ 1 +  Machine + (1|Worker) + (0+Machine|Worker)), machines)
            # # as a cheat for comparing the covariance matrices, we use packages
            # @test first(rlmm.rePCA) ≈ first(jlmm.rePCA) atol=1e-3           
        end
        
        

    end
end
