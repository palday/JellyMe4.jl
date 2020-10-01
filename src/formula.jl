function convert_julia_to_r(f::StatsModels.FormulaTerm)::AbstractString
    rhs = Vector{String}(undef, length(f.rhs))
    
    for (idx, trm) in enumerate(f.rhs)

        if trm isa MixedModels.ZeroCorr
            for tt in trm.term.lhs.terms
                if tt isa CategoricalTerm 
                    size(tt.contrasts.matrix, 2) == 2 || 
                        throw(ArgumentError("zerocorr for factors with more than 2 levels not supported")) 
                end 
            end
        
            strtrm = replace(string(trm), "MixedModels.ZeroCorr((" =>"")
            offset = match(r"\|.*\)\)", strtrm).offset
            grp = strip(strtrm[offset+1:end-2])
            pred = strip(strtrm[begin:offset-1])
            septerms = ("""($(tt == "1" ? tt : "0 + $tt") | $grp)""" for tt in strip.(split(pred, r"[+*]")))
            
            rhs[idx] = foldl((x,y) -> "$x + $y", septerms)    
        else
            rhs[idx] = string(trm)    
        end
        # fulldummy
        rhs[idx] = replace(rhs[idx], r"fulldummy\(([^)]*)\)" => s"dummy(\1,levels(\1))")
    end
    
    
    formula = "$(f.lhs) ~ $(foldl((x,y) -> "$x + $y", string.(rhs)))"
    
    formula = replace(formula, ":(log" => "(log")
    formula = replace(formula, ":(exp" => "(exp")
    # zscore won't work here because the operations within a formula are broadcast
    # and never see the whole array
    # formula = replace(formula, ":(zscore" => "(scale")

    occursin(":(", formula) && throw(ArgumentError("Formula contains a transformation."))

    formula = replace(formula, "&" => ":")
end


