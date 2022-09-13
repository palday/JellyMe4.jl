using CategoricalArrays
import LinearAlgebra: pinv
using DataFrames
using MixedModels: getθ, nranef
import RCall: rcopy,
              reval,
              @R_str

"""
    categorical!(df::DataFrame, cols::Vector{Symbol})
    categorical!(df::DataFrame, col::Symbol)

Convert `cols` of `df` to CategoricalArrays.

This replaces functionality previously part of DataFrames.jl and
is intentionally not exported to avoid type piracy.
"""
function categorical!(df::DataFrame, cols::Vector{Symbol})
    return transform!(df, (cc => categorical for cc in cols)...; renamecols=false)
end
function categorical!(df::DataFrame, col::Symbol)
    return transform!(df, col => categorical; renamecols=false)
end

_guarantee_array(val::AbstractArray) = val
_guarantee_array(val) = [val]

"""
    rcopyarray(robj)

Copy an object from R, guaranteeing that it is an array.

This is useful to undo RCall's automatic conversion of size 1 arrays to scalars.
"""
rcopyarray(robj) = _guarantee_array(rcopy(robj))

"""
get the contrasts from an R dataframe
"""
function get_r_contrasts(rdf)
    data = rcopy(rdf)
    # get categorical columns
    cnames = [c for c in propertynames(data) if typeof(data[!, c]) <: CategoricalArray]
    return Dict(c => HypothesisCoding(pinv(rcopy(R"contrasts($(rdf[c]))"));
                                      labels=rcopyarray(R"colnames(contrasts($(rdf[c])))"))
                for c in cnames)
end

"""
set the contrasts on an R dataframe
this of course has a side effect: it changes the R dataframe
"""

function set_r_contrasts!(rdfname, formula)
    fixefform = first(formula.rhs)

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
        elseif typeof(tt) <: InteractionTerm &&
               any(typeof(tx) <: CategoricalTerm for tx in tt.terms)
            # do nothing -- unless you're doing something really crazy,
            # then this should be handled by the coding of the first-order terms
            # if you are doing something that crazy, then you know enough linear algebra
            # to copy and interpret the relevant matrices directly
            @info "contrasts on interaction terms are assumed to decompose into products of contrasts on the first-order terms" maxlog = 1
            @info "(if you don't know what that means, you're probably fine)" maxlog = 1
        end
    end
    return rdfname
end

#####
##### See GitHub issue #38 on ordering
#####

function _reorder_theta_from_lme4(θlme4, model)
    reterms = model.reterms

    # not the most efficient procedure, but it's not too hard

    # put them in lme4 order
    rr = sort(reterms; rev=true, by=x -> length(x.levels))
    # compute the permutation from lme4 to julia
    reperm = sort(1:length(rr); rev=true, by=x -> nranef(rr[x]))

    # extract the individual theta blocks
    start = 1
    θs = []
    for r in rr
        finish = start + length(getθ(r)) - 1
        push!(θs, θlme4[start:finish])
        start = finish + 1
    end

    # permute
    θjulia = vcat(θs[reperm]...)
    return θjulia
end

function _reorder_theta_to_lme4(model)
    # lme4 sorts the order of RE terms based on nlevels of the grouping var
    reperm = sort(1:length(model.reterms);
                  rev=true,
                  by=x -> length(model.reterms[x].levels))

    return vcat(getθ.(model.reterms)[reperm]...)
end
