include("set_up_tests.jl")

@testset ExtendedTestSet "merMod" include("merMod.jl")
@testset ExtendedTestSet "lmerMod" include("lmerMod.jl")
@testset ExtendedTestSet "glmerMod" include("glmerMod.jl")
