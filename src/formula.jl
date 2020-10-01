function convert_julia_to_r(f::StatsModels.FormulaTerm)::AbstractString
    formula = string(f)

    formula = replace(formula, ":(log" => "(log")
    formula = replace(formula, ":(exp" => "(exp")
    # zscore won't work here because the operations within a formula are broadcast
    # and never see the whole array
    # formula = replace(formula, ":(zscore" => "(scale")

    occursin(":(", formula) && throw(ArgumentError("Formula contains a transformation."))

    occursin("zerocorr", lowercase(formula)) && throw(ArgumentError("zerocorr not supported"))
    formula = replace(formula, "&" => ":")
end
