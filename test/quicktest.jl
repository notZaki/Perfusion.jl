using BenchmarkTools, NumericalIntegration
using Perfusion

t = collect(1:600) ./ 60
cp = aif_georgiou(t .- 1)

test_params = (kt = 0.25, kep = 0.5)
ct = model_tofts(t = t, cp = cp, parameters = test_params)
ref_params = (kt = 0.14, kep = 1.0)
crr = model_tofts(t = t, cp = cp, parameters = ref_params)

rel_kt = test_params.kt / ref_params.kt
rel_ve = (test_params.kt / test_params.kep) / (ref_params.kt / ref_params.kep)
rel_kt = round(rel_kt, digits = 3)
rel_ve = round(rel_ve, digits = 3)

num_noisy_replications = 100
σ = 0.1
noisy_ct = zeros(num_noisy_replications, length(t))
for idx = 1:num_noisy_replications
    noisy_ct[idx, :] .= ct .+ σ .* randn(length(ct))
end

estimates = fit_model(:referenceregion, :lls, t = t, crr = crr, ct = noisy_ct).estimates
std(estimates.rel_kt)
std(estimates.rel_ve)
std(estimates.kep)

estimates = fit_model(:constrained_referenceregion, :lls, t = t, crr = crr, ct = noisy_ct).estimates
mean(estimates.rel_kt)
mean(estimates.rel_ve)
mean(estimates.kep)
