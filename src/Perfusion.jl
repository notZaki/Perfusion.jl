module Perfusion

include("utils.jl")
export @extract, interquartile_mean, percent_error, make_folder, readpath

using LinearAlgebra: norm
using LsqFit
using Statistics: mean
include("relaxation.jl")
export spgr,
       concentration_to_R1,
       concentration_to_signal,
       signal_to_concentration,
       signal_to_R1
export fit_relaxation, fit_relaxation_nls, fit_relaxation_despot, fit_relaxation_novifast

using SpecialFunctions: gamma
include("aif.jl")
export aif_biexponential, aif_fritzhansen, aif_weinmann
export aif_orton1, aif_orton2, aif_orton3
export aif_georgiou, aif_parker

using LsqFit
using NamedTupleTools: select
using NumericalIntegration: cumul_integrate, TrapezoidalFast
using Statistics: mean, quantile, std
include("model_fitting.jl")
export model_tofts, model_exchange, model_filtration, model_uptake
export fit_model
export fit_tofts_nls, fit_tofts_lls
export fit_extendedtofts_nls, fit_extendedtofts_lls
export fit_uptake_nls, fit_uptake_lls
export fit_exchange_nls, fit_exchange_lls
export fit_filtration_nls, fit_filtration_lls
export fit_rrm_nls, fit_rrm_lls, fit_crrm_nls, fit_crrm_lls
export fit_errm_lls, fit_cerrm_lls
export fit_rrift, fit_rrift_with_crrm, fit_rrift_with_cerrm

end # module
