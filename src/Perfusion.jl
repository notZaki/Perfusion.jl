module Perfusion

macro extract(varnames, namedtuple)
    ex = Expr(:block)
    ex.args = [:($(esc(var)) = getindex($(esc(namedtuple)), $(esc(QuoteNode(var))))) for var in varnames.args]
    ex
end

using SpecialFunctions: gamma
include("aif.jl")
export aif_biexponential, aif_fritzhansen, aif_weinmann
export aif_orton1, aif_orton2, aif_orton3
export aif_georgiou, aif_parker

using NamedTupleTools: select
using NumericalIntegration: cumul_integrate, TrapezoidalFast
using LsqFit
include("model_fitting.jl")
export model_tofts, model_exchange, model_filtration, model_uptake
export fit_model
export fit_tofts_nls, fit_tofts_lls
export fit_extendedtofts_nls, fit_extendedtofts_lls
export fit_uptake_nls, fit_uptake_lls
export fit_exchange_nls, fit_exchange_lls
export fit_filtration_nls, fit_filtration_lls

end # module
