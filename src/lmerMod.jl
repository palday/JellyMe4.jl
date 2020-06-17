import MixedModels: LinearMixedModel,
                    setθ!,
                    updateL!
import RCall: rcopy,
              RClass,
              rcopytype,
              reval,
              S4Sxp,
              sexp,
              protect,
              unprotect,
              sexpclass,
              @rput,
              @rget,
              @R_str

# if RCall is available, then so is DataFrames
import DataFrames: DataFrame
import Tables: ColumnTable

# from R
# note that weights are not extracted
# TODO: document weights issue and warn
function rcopy(::Type{LinearMixedModel}, s::Ptr{S4Sxp})
    # this only extracts the name within the call, not the actual weights
    try
        wts = rcopy(s[:call][:weights])
        @error "weights are not supported"
    catch err
        if !isa(err, BoundsError) # something we weren't expecting
            throw(err)
        end
        # no weights defined, we continue on our way
    end
    f = rcopy(s[:call][:formula])
    data = rcopy(s[:frame])
    contrasts = get_r_contasts(s[:frame])

    θ = rcopy(s[:theta])
    reml = rcopy(s[:devcomp][:dims][:REML]) ≠ 0

    m = LinearMixedModel(f, data, contrasts=contrasts)
    m.optsum.REML = reml
    m.optsum.feval = rcopy(s[:optinfo][:feval])
    try
        m.optsum.final = rcopy(s[:optinfo][:val])
    catch err
        if isa(err, MethodError)
            # this happens if θ has length one, i.e. a single scalar RE
            m.optsum.final = [rcopy(s[:optinfo][:val])]
            θ = [θ]
        else
            throw(err)
        end
    end
    m.optsum.optimizer = Symbol("$(rcopy(s[:optinfo][:optimizer])) (lme4)")
    m.optsum.returnvalue = Bool(rcopy(s[:optinfo][:conv][:opt])) ? :FAILURE : :SUCCESS
    m.optsum.fmin = reml ? rcopy(s[:devcomp][:cmp][:REML]) : rcopy(s[:devcomp][:cmp][:dev])
    updateL!(setθ!(m, θ))
end

rcopytype(::Type{RClass{:lmerMod}}, s::Ptr{S4Sxp}) = LinearMixedModel

# TODO: fix some conversions -- Julia->R->Julia roundtrip currently due to
#        ERROR: REvalError: Error in function (x, value, pos = -1, envir = as.environment(pos), inherits = FALSE,  :
#          SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
function sexp(::Type{RClass{:lmerMod}}, x::Tuple{LinearMixedModel{T}, DataFrame}) where T
    m, tbl = x
    if !isempty(m.sqrtwts)
        @error "weights are not currently supported"
    end

    m.optsum.feval > 0 || throw(ArgumentError("Model must be fitted"))

    # should we assume the user is smart enough?
    reval("library(lme4)")

    jellyme4_data = tbl
    formula = convert_julia_to_r(m.formula)

    θ = m.θ
    rsteps = 1
    REML = m.optsum.REML ? "TRUE" : "FALSE"
    jellyme4_theta = m.optsum.final
    fval = m.optsum.fmin
    feval = m.optsum.feval
    conv = m.optsum.returnvalue == :SUCCESS ? 0 : 1
    optimizer = String(m.optsum.optimizer)
    message = "fit with MixedModels.jl"
    # yes, it overwrites any variable named data, but you shouldn't be naming
    # your variables that anyway!
    @rput jellyme4_data
    @rput jellyme4_theta

    r = """
    jellyme4_mod <- lmer(formula = $(formula),
                           data=jellyme4_data,
                           REML=$(REML),
                           control=lmerControl(optimizer="nloptwrap",
                                                optCtrl=list(maxeval=$(rsteps)),
                                    calc.derivs=FALSE),
                            start=list(theta=jellyme4_theta))
     jellyme4_mod@optinfo\$feval <- $(feval)
     jellyme4_mod@optinfo\$message <- "$(message)"
     jellyme4_mod@optinfo\$optimizer <- "$(optimizer)"
     jellyme4_mod
    """
    @debug r
    r = reval(r)
    r = protect(sexp(r))
    unprotect(1)
    r
end

sexpclass(x::Tuple{LinearMixedModel{T}, DataFrame}) where T = RClass{:lmerMod}

# generalize to ColumnTable, which is what MixedModels actually requires
function sexp(ss::Type{RClass{:lmerMod}}, x::Tuple{LinearMixedModel{T}, ColumnTable}) where T
    m, t  = x
    sexp(ss, Tuple([m, DataFrame(t)]))
end

sexpclass(x::Tuple{LinearMixedModel{T}, ColumnTable}) where T = RClass{:lmerMod}
