module Perfusion

# Arterial input functions

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

end # module
