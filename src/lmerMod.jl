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

    # these try blocks should probably be changed to an examination of the indices
    # this only extracts the name within the call, not the actual weights
    try
        wts = rcopy(s[:call][:weights])
        @error "weights are not supported"
    catch err
        if !isa(err, BoundsError) # this is the error we were expecting
            throw(err)
        end
        # no weights defined, we continue on our way
    end

    try
        contrasts = rcopy(s[:call][:contrasts])
        @error "Contrasts must be specified in the dataframe, not the lmer() call"
    catch err
        if !isa(err, BoundsError) # this is the error we were expecting
            throw(err)
        end
        # no extra contrasts defined, we continue on our way
    end

    f = rcopy(s[:call][:formula])
    data = rcopy(s[:frame])
    contrasts = get_r_contrasts(s[:frame])

    θ = rcopyarray(s[:theta])
    reml = rcopy(s[:devcomp][:dims][:REML]) ≠ 0

    m = LinearMixedModel(f, data, contrasts=contrasts)
    m.optsum.REML = reml
    m.optsum.feval = rcopy(s[:optinfo][:feval])
    # I'm wondering if this be filled in from the Julia side
    m.optsum.final = rcopyarray(s[:optinfo][:val])
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

    set_r_contrasts!("jellyme4_data", m.formula)

    r = """
    jellyme4_mod <- lmer(formula = $(formula),
                           data=jellyme4_data,
                           REML=$(REML),
                           control=lmerControl(optimizer="nloptwrap",
                                                optCtrl=list(maxeval=$(rsteps)),
                                    calc.derivs=FALSE,
                                    check.nobs.vs.nRE= "warning"),
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
    sexp(ss, (m, DataFrame(t)))
end

sexpclass(x::Tuple{LinearMixedModel{T}, ColumnTable}) where T = RClass{:lmerMod}
