include("./dce_models/tofts.jl")
include("./dce_models/uptake.jl")
include("./dce_models/exchange.jl")
include("./dce_models/filtration.jl")
include("./dce_models/reference_region.jl")

function fit_model(modelname, fitmethod = :default; kwargs...)
    return model_dict[modelname][fitmethod](; kwargs...)
end

const model_dict = Dict{Symbol,Dict{Symbol,Function}}(
    :tofts => Dict{Symbol,Function}(
        :lls => fit_tofts_lls,
        :nls => fit_tofts_nls,
        :default => fit_tofts_nls,
    ),
    :extendedtofts => Dict{Symbol,Function}(
        :lls => fit_extendedtofts_lls,
        :nls => fit_extendedtofts_nls,
        :default => fit_extendedtofts_nls,
    ),
    :uptake => Dict{Symbol,Function}(
        :lls => fit_uptake_lls,
        :nls => fit_uptake_nls,
        :default => fit_uptake_nls,
    ),
    :exchange => Dict{Symbol,Function}(
        :lls => fit_exchange_lls,
        :nls => fit_exchange_nls,
        :default => fit_exchange_nls,
    ),
    :filtration => Dict{Symbol,Function}(
        :lls => fit_filtration_lls,
        :nls => fit_filtration_nls,
        :default => fit_filtration_nls,
    ),
    :rrm => Dict{Symbol,Function}(
        :lls => fit_rrm_lls,
        :nls => fit_rrm_nls,
        :default => fit_rrm_nls,
    ),
    :crrm => Dict{Symbol,Function}(
        :lls => fit_crrm_lls,
        :nls => fit_crrm_nls,
        :default => fit_crrm_nls,
    ),
    :errm => Dict{Symbol,Function}(
        :lls => fit_errm_lls,
        :default => fit_errm_lls,
    ),
    :cerrm => Dict{Symbol,Function}(
        :lls => fit_cerrm_lls,
        :default => fit_cerrm_lls,
    ),
    :rrift => Dict{Symbol,Function}(
        :crrm => fit_crrm_with_rrift,
        :cerrm => fit_cerrm_with_rrift,
        :default => fit_crrm_with_rrift,
    ),
)
