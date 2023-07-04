import MixedModels: MixedModel
import Tables: ColumnTable
import DataFrames: DataFrame
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

function sexp(::Type{RClass{:merMod}}, x::MixedModel{T}) where {T}
    throw(ArgumentError("You must use a tuple with the data -- (MixedModel, DataFrame) -- not just a model"))
end

sexpclass(x::MixedModel{T}) where {T} = RClass{:merMod}

# we could in theory support the other ordering by re-ordering and passing on
# but that introduces additional maintainance work and encourages messy style
function sexp(::Type{RClass{:merMod}}, x::Tuple{Union{DataFrame,ColumnTable},MixedModel})
    throw(ArgumentError("The order in your Tuple is reversed. It should be (Model, Data)"))
end

sexpclass(x::Tuple{Union{DataFrame,ColumnTable},MixedModel}) = RClass{:merMod}
