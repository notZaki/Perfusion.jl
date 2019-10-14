function expconv(A::AbstractVector, B::Number, t::AbstractVector)
    # Returns f = A ConvolvedWith exp(-B t)
    # Based on Flouri et al. (2016) MRM 76(3), doi: 10.1002/mrm.25991
    @assert length(A) == length(t)
    f = zeros(length(t))
    for i in 2:length(t)
        x = B * (t[i] - t[i-1])
        dA = (A[i] - A[i-1]) / x
        E = exp(-x)
        E0 = 1 - E
        E1 = x - E0
        f[i] = E * f[i-1] + A[i-1] * E0 + dA * E1
    end
    return f ./ B
end

function model_tofts(; t::AbstractVector, parameters::NamedTuple, cp::AbstractVector)
    @extract (kt, kep) parameters
    ct = kt * expconv(cp, kep, t)
    if haskey(parameters, :vp)
        ct .+= parameters.vp * cp
    end
    return ct
end

function model_exchange(; t::AbstractVector, parameters::NamedTuple, ca::AbstractVector)
    if all(haskey(parameters, key) for key in (:Tp, :Te))
        @extract (ve, vp) parameters
        fp = vp / parameters.Tp
        ps = ve / parameters.Te
    else
        @extract (fp, ps, ve, vp) parameters
    end
    Tp = vp / fp
    Te = ve / ps
    T = (vp + ve) / fp
    Tplus = (T + Te + sqrt((T + Te)^2 - 4 * Tp * Te)) / 2
    Tminus = (T + Te - sqrt((T + Te)^2 - 4 * Tp * Te)) / 2
    ct = ((T - Tminus) / (Tplus - Tminus)) .* expconv(ca, 1/Tplus, t)
    ct .+= ((Tplus - T) / (Tplus - Tminus)) .* expconv(ca, 1/Tminus, t)
    ct .*= fp
    return ct
end

function fit_exchange_lls(; t::AbstractVector, ca::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    α, β, γ, fp = (fill(NaN, volume_size...) for _=1:4)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 4)
    M[:,4] .= cumul_integrate(t, ca, TrapezoidalFast())
    M[:,3] .= cumul_integrate(t, M[:,4], TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:,2] .= -cumul_integrate(t, ct[idx,:], TrapezoidalFast())
        M[:,1] .= cumul_integrate(t, M[:,2], TrapezoidalFast())
        α[idx], β[idx], γ[idx], fp[idx] = M \ ct[idx,:]
    end
    T  = @. γ / (α * fp)
    Te = @. (β / α) - T
    Tp = @. 1 / (α * Te)
    vp = @. fp * Tp
    ve = @. fp * (T - Tp)
    ps = @. ve / Te
    return(estimates=(fp=fp, ps=ps, ve=ve, vp=vp, T=T, Te=Te, Tp=Tp), dummy=0)
end

function fit_exchange_nls(; t::AbstractVector, ca::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    fp, ps, vp, ve = (fill(NaN, volume_size...) for _=1:4)
    resolved_mask = resolve_mask_size(mask, volume_size)
    model(x, p) = model_exchange(t=x, ca=ca, parameters=(fp=p[1], ps=p[2], ve=p[3], vp=p[4]))
    lls_estimates = fit_exchange_lls(t=t, ca=ca, ct=ct).estimates
    init_fp, init_ps, init_ve, init_vp = select(lls_estimates, (:fp, :ps, :ve, :vp))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalue = [init_fp[idx], init_ps[idx], init_ve[idx], init_vp[idx]]
        fp[idx], ps[idx], ve[idx], vp[idx] = curve_fit(model, t, ct[idx, :], initialvalue).param
    end
    T  = @. (vp + ve) / fp
    Tp = @. vp / fp
    Te = @. ve / ps
    return(estimates=(fp=fp, ps=ps, ve=ve, vp=vp, T=T, Tp=Tp, Te=Te), dummy=0)
end

function model_filtration(; t::AbstractVector, parameters::NamedTuple, ca::AbstractVector)
    if all(haskey(parameters, key) for key in (:Tp, :Te))
        @extract (ve, vp) parameters
        fp = vp / parameters.Tp
        ps = ve / parameters.Te
    else
        @extract (fp, ps, ve, vp) parameters
    end
    Tminus = vp/fp
    Tplus = ve/ps
    T = (vp+ve)/fp
    ct = ((T - Tminus) / (Tplus - Tminus)) .* expconv(ca, 1/Tplus, t)
    ct .+= ((Tplus - T) / (Tplus - Tminus)) .* expconv(ca, 1/Tminus, t)
    ct .*= fp
    return ct
end

function fit_filtration_lls(; t::AbstractVector, ca::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    α, β, γ, fp = (fill(NaN, volume_size...) for _=1:4)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 4)
    M[:,4] .= cumul_integrate(t, ca, TrapezoidalFast())
    M[:,3] .= cumul_integrate(t, M[:,4], TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:,2] .= -cumul_integrate(t, ct[idx,:], TrapezoidalFast())
        M[:,1] .= cumul_integrate(t, M[:,2], TrapezoidalFast())
        α[idx], β[idx], γ[idx], fp[idx] = M \ ct[idx,:]
    end
    root_term = @. β^2 - 4 * α
    @. root_term[root_term < 0] = 0
    Tp = @. (β - sqrt(root_term)) / (2 * α)
    Te = @. (β + sqrt(root_term)) / (2 * α)
    T  = @. γ / (α * fp)
    vp = @. fp * Tp
    ve = @. fp * (T - Tp)
    ps = @. ve / Te
    return(estimates=(fp=fp, ps=ps, ve=ve, vp=vp, T=T, Te=Te, Tp=Tp), dummy=0)
end

function fit_filtration_nls(; t::AbstractVector, ca::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    fp, ps, vp, ve = (fill(NaN, volume_size...) for _=1:4)
    resolved_mask = resolve_mask_size(mask, volume_size)
    model(x, p) = model_filtration(t=x, ca=ca, parameters=(fp=p[1], ps=p[2], ve=p[3], vp=p[4]))
    lls_estimates = fit_filtration_lls(t=t, ca=ca, ct=ct).estimates
    init_fp, init_ps, init_ve, init_vp = select(lls_estimates, (:fp, :ps, :ve, :vp))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalue = [init_fp[idx], init_ps[idx], init_ve[idx], init_vp[idx]]
        fp[idx], ps[idx], ve[idx], vp[idx] = curve_fit(model, t, ct[idx, :], initialvalue).param
    end
    T  = @. (vp + ve) / fp
    Tp = @. vp / fp
    Te = @. ve / ps
    @. T[fp == 0] = 0
    @. Tp[fp == 0] = 0
    @. Te[ps == 0] = 0
    return(estimates=(fp=fp, ps=ps, ve=ve, vp=vp, T=T, Tp=Tp, Te=Te), dummy=0)
end

function model_uptake(; t::AbstractVector, parameters::NamedTuple, ca::AbstractVector)
    @extract (fp, ps, vp) parameters
    Tp = vp / (fp + ps)
    E = ps / (fp + ps)
    ca_conv_1 = zeros(length(t))
    for idx in 2:length(ca)
        ca_conv_1[idx] = ca_conv_1[idx-1] + (t[idx]-t[idx-1]) * ca[idx]
    end
    ct = fp .* ((1-E) .* expconv(ca, 1/Tp, t) .+ E .* ca_conv_1)
    return ct
end

function fit_uptake_lls(; t::AbstractVector, ca::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    x1, fp, x3, ps, vp = (fill(NaN, volume_size...) for _=1:5)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 3)
    M[:,2] = cumul_integrate(t, ca, TrapezoidalFast())
    M[:,3] = cumul_integrate(t, M[:,2], TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 1] = -cumul_integrate(t, ct[idx,:], TrapezoidalFast())
        (x1[idx], fp[idx], x3[idx]) = M \ ct[idx,:]
    end
    denum = @. x1 * fp - x3
    @. ps = fp * x3 / denum
    @. vp = fp^2 / denum
    @. ps[denum == 0] = 0
    @. vp[denum == 0] = 0
    return(estimates=(fp=fp, ps=ps, vp=vp), dummy=0)
end

function fit_uptake_nls(; t::AbstractVector, ca::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    fp, ps, vp = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)
    model(x, p) = model_uptake(t=x, ca=ca, parameters=(fp=p[1], ps=p[2], vp=p[3]))
    lls_estimates = fit_uptake_lls(t=t, ca=ca, ct=ct).estimates
    init_fp, init_ps, init_vp = select(lls_estimates, (:fp, :ps, :vp))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalue = [init_fp[idx], init_ps[idx], init_vp[idx]]
        (fp[idx], ps[idx], vp[idx]) = curve_fit(model, t, ct[idx, :], initialvalue).param
    end
    return(estimates=(fp=fp, ps=ps, vp=vp), dummy=0)
end

function fit_tofts_nls(; t::AbstractVector, cp::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(cp) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    kt, kep = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)

    model(x, p) = model_tofts(t=x, cp=cp, parameters=(kt=p[1], kep=p[2]))
    initialvalues = [0.01, 0.01]
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        (kt[idx], kep[idx]) = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    return(estimates=(kt=kt, kep=kep), dummy=0)
end

function fit_tofts_lls(; t::AbstractVector, cp::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(cp) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    kt, kep = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 2)
    M[:,1] = cumul_integrate(t, cp, TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, ct[idx,:], TrapezoidalFast())
        (kt[idx], kep[idx]) = M \ ct[idx,:]
    end
    return(estimates=(kt=kt, kep=kep), dummy=0)
end

function fit_extendedtofts_lls(;t::AbstractVector, cp::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(cp) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    kt, kep, vp = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 3)
    M[:,1] = cumul_integrate(t, cp, TrapezoidalFast())
    M[:,3] = cp
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, ct[idx,:], TrapezoidalFast())
        (kt[idx], kep[idx], vp[idx]) = M \ ct[idx,:]
    end
    # Apply correction because fit actually returns: kt + kep * vp
    @. kt = kt - (kep * vp)
    return(estimates=(kt=kt, kep=kep, vp=vp), dummy=0)
end

function fit_extendedtofts_nls(; t::AbstractVector, cp::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(cp) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    kt, kep, vp = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)

    model(x, p) = model_tofts(t=x, cp=cp, parameters=(kt=p[1], kep=p[2], vp=p[3]))
    initialvalues = [0.01, 0.01, 0.01]
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        (kt[idx], kep[idx], vp[idx]) = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    return(estimates=(kt=kt, kep=kep, vp=vp), dummy=0)
end

function model_referenceregion(; t::AbstractVector, parameters::NamedTuple, crr::AbstractVector)
    @extract (rel_kt, kep, kep_rr) parameters
    ct = rel_kt .* crr .+ (rel_kt * (kep_rr - kep)) .* expconv(crr, kep, t)
    return ct
end

function fit_referenceregion_nls(; t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(crr) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    rel_kt, kep, kep_rr = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)
    lls_estimates = fit_referenceregion_lls(t=t, crr=crr, ct=ct).estimates
    init_rel_kt, init_kep, init_kep_rr = select(lls_estimates, (:rel_kt, :kep, :kep_rr))
    model(x, p) = model_referenceregion(t=x, crr=crr,
        parameters=(rel_kt=p[1], kep=p[2], kep_rr=p[3]))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalues = [init_rel_kt[idx], init_kep[idx], init_kep_rr[idx]]
        rel_kt[idx], kep[idx], kep_rr[idx] = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return(estimates=(rel_kt=rel_kt, rel_ve=rel_ve, kep=kep, kep_rr=kep_rr), dummy=0)
end

function fit_referenceregion_lls(;t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true)
    @assert length(t) == length(crr) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    x1, x2, x3 = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 3)
    M[:,1] = crr
    M[:,2] = cumul_integrate(t, crr, TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 3] = -cumul_integrate(t, ct[idx,:], TrapezoidalFast())
        x1[idx], x2[idx], x3[idx] = M \ ct[idx,:]
    end
    # Apply correction because fit actually returns: [rel_kt, kt/ve_rr, kep]
    rel_kt = x1
    rel_ve = x2 ./ x3
    kep = x3
    kep_rr = x2 ./ x1
    return(estimates=(rel_kt=rel_kt, rel_ve=rel_ve, kep=kep, kep_rr=kep_rr), dummy=0)
end

function fit_constrained_referenceregion_nls(; t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true, kep_rr=0.0)
    @assert length(t) == length(crr) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    rel_kt, kep = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)
    if kep_rr <= 0
        rrm_estimates = fit_referenceregion_nls(t=t, crr=crr, ct=ct, mask=mask).estimates
        viable_estimates = positive_only_mask(rrm_estimates)
        kep_rr = interquartile_mean(rrm_estimates.kep_rr[viable_estimates])
    end
    lls_estimates = fit_constrained_referenceregion_lls(t=t, crr=crr, ct=ct, kep_rr=kep_rr).estimates
    init_rel_kt, init_kep = select(lls_estimates, (:rel_kt, :kep))
    model(x, p) = model_referenceregion(t=x, crr=crr,
        parameters=(rel_kt=p[1], kep=p[2], kep_rr=kep_rr))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalues = [init_rel_kt[idx], init_kep[idx]]
        rel_kt[idx], kep[idx] = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return(estimates=(rel_kt=rel_kt, rel_ve=rel_ve, kep=kep, kep_rr=kep_rr), dummy=0)
end

function fit_constrained_referenceregion_lls(; t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true, kep_rr=0.0)
    @assert length(t) == length(crr) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    rel_kt, kep = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)

    if kep_rr <= 0
        rrm_estimates = fit_referenceregion_lls(t=t, crr=crr, ct=ct, mask=mask).estimates
        viable_estimates = positive_only_mask(rrm_estimates)
        kep_rr = interquartile_mean(rrm_estimates.kep_rr[viable_estimates])
    end

    M = zeros(num_timepoints, 2)
    M[:,1] = crr + kep_rr * cumul_integrate(t, crr, TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, ct[idx,:], TrapezoidalFast())
        rel_kt[idx], kep[idx] = M \ ct[idx,:]
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return(estimates=(rel_kt=rel_kt, rel_ve=rel_ve, kep=kep, kep_rr=kep_rr), dummy=0)
end

function positive_only_mask(x::NamedTuple)
    valid_keys = keys(x)
    mask = trues(size(x[valid_keys[1]]))
    for array in x
        @. mask = mask & (array > 0)
    end
    return mask
end

function interquartile_mean(x::AbstractVector)
    if std(x)>1e-3
        quartiles = quantile(x, [0.25, 0.75])
        # `<=` is safer than `<` because latter can creat empty array if x is short
        interquartile_x = @. x[quartiles[1] <= x <= quartiles[2]]
    else
        interquartile_x = x
    end
    return mean(interquartile_x)
end

function resolve_mask_size(mask, desired_size)
    if size(mask) == desired_size
        return mask .> 0
    elseif size(mask) == ()
        return repeat([mask .> 0], desired_size...)
    else
        error("Mask size: $(size(mask)) does not match input size $(desired_size)")
    end
end

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
