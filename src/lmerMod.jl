# from R
# note that weights are not extracted
# TODO: document weights issue and warn
function RCall.rcopy(::Type{LinearMixedModel}, s::Ptr{S4Sxp})

    # these try blocks should probably be changed to an examination of the indices
    # this only extracts the name within the call, not the actual weights
    try
        wts = rcopy(s[:call][:weights])
        @error "weights are not supported"
    catch err
        if !isa(err, BoundsError) # this is the error we were expecting
            rethrow(err)
        end
        # no weights defined, we continue on our way
    end

    try
        contrasts = rcopy(s[:call][:contrasts])
        @error "Contrasts must be specified in the dataframe, not the lmer() call"
    catch err
        if !isa(err, BoundsError) # this is the error we were expecting
            rethrow(err)
        end
        # no extra contrasts defined, we continue on our way
    end

    # for some reason this doesn't always give a formula with lmerTest
    #f = rcopy(s[:call][:formula])
    f = rcopy(R"as.formula($(s)@call$formula)")
    data = rcopy(s[:frame])
    contrasts = get_r_contrasts(s[:frame])
    θ = rcopyarray(s[:theta])
    reml = rcopy(s[:devcomp][:dims][:REML]) ≠ 0

    m = LinearMixedModel(f, data; contrasts=contrasts)

    if length(θ) != length(m.θ)
        @error """You're probably using || in R with a categorical variable,
                  whose translation is currently unsupported with recent MixedModels releases"""
        throw(ArgumentError("Parameter vectors in R and Julia are different sizes."))
    end

    θ = _reorder_theta_from_lme4(θ, m)

    m.optsum.REML = reml
    m.optsum.feval = rcopy(s[:optinfo][:feval])
    # I'm wondering if this be filled in from the Julia side
    m.optsum.final = rcopyarray(s[:optinfo][:val])
    m.optsum.optimizer = Symbol("$(rcopy(s[:optinfo][:optimizer])) (lme4)")
    m.optsum.returnvalue = rcopy(s[:optinfo][:conv][:opt]) == 0 ? :FAILURE : :SUCCESS
    m.optsum.fmin = reml ? rcopy(s[:devcomp][:cmp][:REML]) : rcopy(s[:devcomp][:cmp][:dev])
    return updateL!(setθ!(m, θ))
end

RCall.rcopytype(::Type{RClass{:lmerMod}}, s::Ptr{S4Sxp}) = LinearMixedModel

# add lmerTest::lmer and afex::lmer_alt support
RCall.rcopytype(::Type{RClass{:lmerModLmerTest}}, s::Ptr{S4Sxp}) = LinearMixedModel

# TODO: fix some conversions -- Julia->R->Julia roundtrip currently due to
#        ERROR: REvalError: Error in function (x, value, pos = -1, envir = as.environment(pos), inherits = FALSE,  :
#          SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
function RCall.sexp(::Type{RClass{:lmerMod}},
                    x::Tuple{LinearMixedModel{T},DataFrame}) where {T}
    m, tbl = x
    if !isempty(m.sqrtwts)
        @error "weights are not currently supported"
    end

    m.optsum.feval > 0 || throw(ArgumentError("Model must be fitted"))

    jellyme4_data = tbl
    formula = convert_julia_to_r(m.formula)

    θ = m.θ
    REML = m.optsum.REML ? "TRUE" : "FALSE"
    jellyme4_theta = _reorder_theta_to_lme4(m)
    fval = m.optsum.fmin
    feval = m.optsum.feval
    # TODO support things besides NLopt
    conv = m.optsum.returnvalue == :SUCCESS ? 0 : 1
    optimizer = String(m.optsum.optimizer)
    message = "fit with MixedModels.jl"
    @rput jellyme4_data
    @rput jellyme4_theta

    set_r_contrasts!("jellyme4_data", m.formula)

    r = """
    jellyme4_mod <- $LMER(formula = $(formula),
                           data=jellyme4_data,
                           REML=$(REML),
                           control=lme4::lmerControl(optimizer=NULL,
                                                     $(join(MERCONTROL_OPTIONS,","))),
                           start=list(theta=jellyme4_theta))
     jellyme4_mod@optinfo\$feval <- $(feval)
     jellyme4_mod@optinfo\$message <- "$(message)"
     jellyme4_mod@optinfo\$optimizer <- "$(optimizer)"
     jellyme4_mod
    """
    r = reval(r)
    r = protect(sexp(r))
    unprotect(1)
    return r
end

function RCall.sexpclass(x::Tuple{LinearMixedModel{T},DataFrame}) where {T}
    return LMER in ("afex::lmer_alt", "lmer_alt") ? RClass{:lmerModLmerTest} :
           RClass{:lmerMod}
end
function RCall.sexp(::Type{RClass{:lmerModLmerTest}},
                    x::Tuple{LinearMixedModel{T},DataFrame}) where {T}
    return sexp(RClass{:lmerMod}, x)
end

# generalize to ColumnTable, which is what MixedModels actually requires
function RCall.sexp(ss::Type{RClass{:lmerMod}},
                    x::Tuple{LinearMixedModel{T},ColumnTable}) where {T}
    m, t = x
    return sexp(ss, (m, DataFrame(t)))
end

RCall.sexpclass(x::Tuple{LinearMixedModel{T},ColumnTable}) where {T} = RClass{:lmerMod}
