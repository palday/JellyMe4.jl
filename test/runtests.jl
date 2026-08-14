include("set_up_tests.jl")

using Aqua

@testset "Aqua" begin
    Aqua.test_all(JellyMe4;
                  ambiguities=false,
                  deps_compat=(check_extras=true,),
                  piracies=(treat_as_own=[MixedModel, LinearMixedModel,
                                         GeneralizedLinearMixedModel],))
end

@testset "utilities" begin
    df = DataFrame(a=["x", "y", "z"], b=["p", "q", "r"])
    categorical!(df, [:a, :b])
    @test df.a isa CategoricalArray
    @test df.b isa CategoricalArray
end

@testset ExtendedTestSet "merMod" include("merMod.jl")
@testset ExtendedTestSet "lmerMod" include("lmerMod.jl")
@testset ExtendedTestSet "glmerMod" include("glmerMod.jl")
