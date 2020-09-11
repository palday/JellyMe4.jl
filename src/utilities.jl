import CategoricalArrays: CategoricalArray
import LinearAlgebra: pinv
using DataFrames
import RCall: rcopy,
              reval,
              @R_str

"""
    get the contrasts from an R dataframe
"""
function get_r_contrasts(rdf)
    data = rcopy(rdf)
    # get categorical columns
    cnames = [c for c in propertynames(data)  if typeof(data[!, c]) <: CategoricalArray]
    Dict(c => HypothesisCoding(pinv(rcopy(R"contrasts($(rdf[c]))")),
                               labels=rcopy(R"colnames(contrasts($(rdf[c])))")) for c in cnames)
end

"""
    set the contrasts on an R dataframe
    this of course has a side effect: it changes the R dataframe
"""

function set_r_contrasts!(rdfname, formula)
    fixefform = first(formula.rhs)

    warned_on_contrasts = false

    for tt in fixefform.terms
        if typeof(tt) <: CategoricalTerm
            R"""
            jellyme_contrasts <- $(tt.contrasts.matrix)
            colnames(jellyme_contrasts) <- $(tt.contrasts.termnames)
            rownames(jellyme_contrasts) <- $(tt.contrasts.levels)
            """
            # we need rdfname to be interpolated normally before the R macros stuff is called
            reval("""
            # we need to be double-plus sure these are factors
            $(rdfname)[, "$(tt.sym)"] <- factor($(rdfname)[, "$(tt.sym)"])
            contrasts($(rdfname)[, "$(tt.sym)"]) <- jellyme_contrasts
            """)
        elseif typeof(tt) <: InteractionTerm && any( typeof(tx) <: CategoricalTerm for tx in tt.terms)
            # do nothing -- unless you're doing something really crazy,
            # then this should be handled by the coding of the first-order terms
            # if you are doing something that crazy, then you know enough linear algebra
            # to copy and interpret the relevant matrices directly
            if !warned_on_contrasts
                @info "contrasts on interaction terms are assumed to decompose into products of contrasts on the first-order terms"
                @info "(if you don't know what that means, you're probably fine)"
                warned_on_contrasts = true
            end
        end
    end
    rdfname
end
