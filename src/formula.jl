import RCall: rcopy,
              RObject,
              LangSxp,
              @rput,
              @rget,
              @R_str

function convert_julia_to_r(f::StatsModels.FormulaTerm)::AbstractString
    rhs = Vector{String}(undef, length(f.rhs))
    
    for (idx, trm) in enumerate(f.rhs)

        if trm isa MixedModels.ZeroCorr
            grp = trm.term.rhs.sym
            terms = Vector{String}(undef, length(trm.term.lhs.terms))
            
            
            for (idx,tt) in enumerate(trm.term.lhs.terms)
                if tt isa CategoricalTerm 
                    size(tt.contrasts.matrix, 2) == 2 || 
                        throw(ArgumentError("zerocorr for factors with more than 2 levels not supported")) 
                end
                
                if tt isa InterceptTerm
                    terms[idx] = "( 1 | $grp)"
                elseif tt.contrasts isa StatsModels.ContrastsMatrix{StatsModels.FullDummyCoding}
                    terms[idx] = "( 0 + dummy($(tt.sym), levels($(tt.sym))) | $grp)"
                else
                    terms[idx] = "( 0+ $(tt.sym) | $grp)"
                end
                     
            end
        
            # strtrm = replace(string(trm), "MixedModels.ZeroCorr((" =>"")
            # offset = match(r"\|.*\)\)", strtrm).offset
            # grp = strip(strtrm[offset+1:end-2])
            # pred = strip(strtrm[begin:offset-1])
            # septerms = ("""($(tt == "1" ? tt : "0 + $tt") | $grp)""" for tt in strip.(split(pred, r"[+*]")))
            
            rhs[idx] = foldl((x,y) -> "$x + $y", terms)    
        elseif trm isa MixedModels.RandomEffectsTerm
            grp = trm.rhs.sym
            terms = Vector{String}(undef, length(trm.lhs.terms))
            
            for (idx,tt) in enumerate(trm.lhs.terms)   
                if tt isa CategoricalTerm && 
                    tt.contrasts isa StatsModels.ContrastsMatrix{StatsModels.FullDummyCoding}
                    terms[idx] = "dummy($(tt.sym), levels($(tt.sym)))"
                else
                    terms[idx] = string(tt)
                end
                     
            end
        
            rhs[idx] = "( " * foldl((x,y) -> "$x + $y", terms) * " | $grp)"
        else
            rhs[idx] = string(trm)    
        end
        # fulldummy
        rhs[idx] = replace(rhs[idx], r"fulldummy\(([^)]*)\)" => s"dummy(\1,levels(\1))")
    end
    
    
    formula = "$(f.lhs) ~ $(foldl((x,y) -> "$x + $y", string.(rhs)))"
    
    formula = replace(formula, ":(log" => "(log")
    formula = replace(formula, ":(exp" => "(exp")
    # should we consider supporting division on the lhs?
    # what about nesting syntax on the rhs?
    # zscore won't work here because the operations within a formula are broadcast
    # and never see the whole array
    # formula = replace(formula, ":(zscore" => "(scale")

    occursin(":(", formula) && throw(ArgumentError("Formula contains a transformation."))

    formula = replace(formula, "&" => ":")
end

# NOTE: this is known to yet to work properly!
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
    # @info vv
    # @info typeof(vv)
    return vv
end

convert_r_to_julia(f::Ptr{LangSxp})::StatsModels.FormulaTerm =
                convert_r_to_julia(RObject{LangSxp}(f))
