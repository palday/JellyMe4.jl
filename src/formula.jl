import RCall: rcopy,
              RObject,
              LangSxp,
              @rput,
              @rget,
              @R_str
              
function convert_julia_to_r(f::StatsModels.FormulaTerm)::AbstractString
    formula = string(f)

    formula = replace(formula, ":(log" => "(log")
    formula = replace(formula, ":(exp" => "(exp")
    # should we consider supporting division on the lhs?
    # what about nesting syntax on the rhs?
    # zscore won't work here because the operations within a formula are broadcast
    # and never see the whole array
    # formula = replace(formula, ":(zscore" => "(scale")

    occursin(":(", formula) && throw(ArgumentError("Formula contains a transformation."))

    occursin("zerocorr", formula) && throw(ArgumentError("zerocorr not supported"))
    formula = replace(formula, "&" => ":")
end

function convert_r_to_julia(f::RObject{LangSxp})::StatsModels.FormulaTerm
    _, lhs, rhs = rcopy(R"as.character($f)")
    cbind =  match(r"cbind\((.*),(.*)\)", lhs);
    if isnothing(cbind)
        # we do weird things because of the way / is translated
        vv = rcopy(Expr, f)
    else
        # we don't need to compute the weights from the cbind formula because we can
        # get them directly from the glmerMod object
        numerator = successes = first(cbind.captures);
        failures = last(cbind.captures);
        denominator = "( $successes + $failures )";
        lhs = "$(numerator) / $(denominator)";
        form = "$lhs ~ $rhs"
        # we do weird things because of the way / is translated
        vv = rcopy(Expr, R"""as.formula($form)""")

    end
    vv = eval(:@formula($vv))
    @info vv
    @info typeof(vv)
    return vv
end

convert_r_to_julia(f::Ptr{LangSxp})::StatsModels.FormulaTerm =
                convert_r_to_julia(RObject{LangSxp}(f))
