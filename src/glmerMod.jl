import MixedModels: GeneralizedLinearMixedModel,
                    setθ!,
                    updateL!
import Distributions: Distribution
import GLM: Link
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
              @rget
# if RCall is available, then so is DataFrames
import DataFrames: DataFrame, categorical!
import CategoricalArrays: CategoricalArray
import Tables: ColumnTable
# # from R
# # note that weights are not extracted
# # TODO: document weights issue and warn
 function rcopy(::Type{GeneralizedLinearMixedModel}, s::Ptr{S4Sxp})
     throw(ArgumentError("Sorry, getting GLMMs from R has not been implemented yet."))
#     # this only extracts the name within the call, not the actual weights
#     try
#         wts = rcopy(s[:call][:weights])
#         @error "weights are not supported"
#     catch err
#         if !isa(err, BoundsError) # something we weren't expecting
#             throw(err)
#         end
#         # no weights defined, we continue on our way
#     end
#     f = rcopy(s[:call][:formula])
#     data = rcopy(s[:frame])
#     θ = rcopy(s[:theta])
#     reml = rcopy(s[:devcomp][:dims][:REML]) ≠ 0
#
#     m = LinearMixedModel(f,data)
#     m.optsum.REML = reml
#     m.optsum.feval = rcopy(s[:optinfo][:feval])
#     try
#         m.optsum.final = rcopy(s[:optinfo][:val])
#     catch err
#         if isa(err, MethodError)
#             # this happens if θ has length one, i.e. a single scalar RE
#             m.optsum.final = [rcopy(s[:optinfo][:val])]
#             θ = [θ]
#         else
#             throw(err)
#         end
#     end
#     m.optsum.optimizer = Symbol("$(rcopy(s[:optinfo][:optimizer])) (lme4)")
#     m.optsum.returnvalue = Bool(rcopy(s[:optinfo][:conv][:opt])) ? :FAILURE : :SUCCESS
#     m.optsum.fmin = reml ? rcopy(s[:devcomp][:cmp][:REML]) : rcopy(s[:devcomp][:cmp][:dev])
#     updateL!(setθ!(m, θ))
end

rcopytype(::Type{RClass{:glmerMod}}, s::Ptr{S4Sxp}) = GeneralizedLinearMixedModel

# from Julia
function sexp(::Type{RClass{:glmerMod}}, x::Tuple{GeneralizedLinearMixedModel{T}, DataFrame}) where T
    m, tbl = x
    m.optsum.feval > 0 || throw(ArgumentError("Model must be fitted"))

    distribution =  Distribution(m.resp)
    family = urfamily = replace(string(distribution), Regex("{.*}") => "")
    # ["Bernoulli","Binomial","Gamma", "Poisson", "InverseGaussian"]
    # R families: binomial, gaussian, Gamma, inverse.gaussian, poisson
    #             quasi, quasibinomial, quasipoisson
    # should we @warn for the way GLMM deviance is defined in lme4?
    if family == "Bernoulli"
        family = "binomial"
    else
        throw(ArgumentError("Family $family is not supported"))
    end

    urlink = string(Link(m.resp))
    link = lowercase(replace(urlink, "Link()" => ""))

    link in ["logit", "probit", "cauchit",
             "log", "identity", "inverse", "sqrt",
             "cloglog"] || throw(ArgumentError("Link $urlink not supported"))

    if length(m.optsum.initial) - length(m.β) < 0
        nAGQ = 0
    else
        nAGQ = m.optsum.nAGQ
    end

    # if !isempty(m.sqrtwts)
    #     @error "weights are not currently supported"
    # end


    # should we assume the user is smart enough?
    reval("library(lme4)")

    jellyme4_data = tbl
    if urfamily == "Bernoulli"
        lhs = m.formula.lhs
        # should we check for PooledArray? brings in another dependency...
        if !(jellyme4_data[!,lhs.sym] isa CategoricalArray)
            @warn "Response was not categorical, converting in place"
            categorical!(jellyme4_data, [lhs.sym])
        end
    end
    formula = convert_julia_to_r(m.formula)

    jellyme4_theta = m.θ
    jellyme4_beta = m.β
    fval = m.optsum.fmin
    feval = m.optsum.feval
    conv = m.optsum.returnvalue == :SUCCESS ? 0 : 1
    optimizer = String(m.optsum.optimizer)
    message = "fit with MixedModels.jl"
    rsteps = 1;

    betastart = nAGQ > 0 ? "fixef=jellyme4_beta, " : ""

    @rput jellyme4_data
    @rput jellyme4_theta
    @rput jellyme4_beta

    set_r_contrasts!("jellyme4_data", m.formula)

    r = """
        jellyme4_mod <- glmer(formula = $(formula),
                               data=jellyme4_data,
                               family=$family(link="$link"),
                               nAGQ=$(nAGQ),
                               control=glmerControl(optimizer="nloptwrap",
                                                    optCtrl=list(maxeval=$(rsteps)),
                                        calc.derivs=FALSE),
                                start=list($(betastart)theta=jellyme4_theta))
         jellyme4_mod@optinfo\$feval <- $(feval)
         jellyme4_mod@optinfo\$message <- "$(message)"
         jellyme4_mod@optinfo\$optimizer <- "$(optimizer)"
         jellyme4_mod
    """
    @debug r
    @warn """Some accuracy is lost in translation for GLMMs. This is fine with plotting, but proceed with caution for inferences.
             You can try letting the model "reconverge" in lme4 with
                 update(model, control=glmerControl(optimizer="nloptwrap", calc.derivs=FALSE))."""
    @info "lme4 handles deviance for GLMMs differently than MixedModels.jl"
    @info "for the correct comparison, examine -2*logLik(RModel) and deviance(JuliaModel)"
    r = reval(r)
    r = protect(sexp(r))
    unprotect(1)
    r
end

sexpclass(x::Tuple{GeneralizedLinearMixedModel{T}, DataFrame}) where T = RClass{:glmerMod}

# generalize to ColumnTable, which is what MixedModels actually requires
function sexp(ss::Type{RClass{:glmerMod}}, x::Tuple{GeneralizedLinearMixedModel{T}, ColumnTable}) where T
    m, t  = x
    sexp(ss, Tuple([m, DataFrame(t)]))
end

sexpclass(x::Tuple{GeneralizedLinearMixedModel{T}, ColumnTable}) where T = RClass{:glmerMod}
