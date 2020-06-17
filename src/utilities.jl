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

function set_r_contrasts(rdf, formula)
    fixefform = first(formula.rhs)

    #cc.matrix
    #cc.termnames
end
