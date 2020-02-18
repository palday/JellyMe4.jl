import MixedModels: MixedModel

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

function sexp(::Type{RClass{:merMod}}, x::MixedModel{T}) where T
    throw(ArgumentError("You must use a Tuple(MixedModel, DataFrame), not just a model"))
end

sexpclass(x::MixedModel{T}) where T = RClass{:merMod}
