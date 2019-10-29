include("./dce_models/tofts.jl")
include("./dce_models/uptake.jl")
include("./dce_models/exchange.jl")
include("./dce_models/filtration.jl")
include("./dce_models/reference_region.jl")

function fit_model(modelname, fitmethod=:nls; kwargs...)
    return model_dict[modelname][fitmethod](; kwargs...)
end

const model_dict = Dict{Symbol, Dict{Symbol, Function}}(
    :tofts => Dict{Symbol, Function}(
        :lls => fit_tofts_lls,
        :nls => fit_tofts_nls
    ),
    :extendedtofts => Dict{Symbol, Function}(
        :lls => fit_extendedtofts_lls,
        :nls => fit_extendedtofts_nls
    ),
    :uptake => Dict{Symbol, Function}(
        :lls => fit_uptake_lls,
        :nls => fit_uptake_nls
    ),
    :exchange => Dict{Symbol, Function}(
        :lls => fit_exchange_lls,
        :nls => fit_exchange_nls
    ),
    :filtration => Dict{Symbol, Function}(
        :lls => fit_filtration_lls,
        :nls => fit_filtration_nls
    ),
    :referenceregion => Dict{Symbol, Function}(
        :lls => fit_referenceregion_lls,
        :nls => fit_referenceregion_nls
    ),
    :constrained_referenceregion => Dict{Symbol, Function}(
        :lls => fit_constrained_referenceregion_lls,
        :nls => fit_constrained_referenceregion_nls
    )
)
