function model_referenceregion(conv::Function=expconv; t::AbstractVector, parameters::NamedTuple, crr::AbstractVector)
    @extract (rel_kt, kep, kep_rr) parameters
    ct = rel_kt .* crr .+ (rel_kt * (kep_rr - kep)) .* conv(crr, kep, t)
    return ct
end

function fit_referenceregion_nls(conv::Function=expconv;
    t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true)

    @assert length(t) == length(crr) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    rel_kt, kep, kep_rr = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)
    lls_estimates = fit_referenceregion_lls(t=t, crr=crr, ct=ct, mask=mask).estimates
    init_rel_kt, init_kep, init_kep_rr = select(lls_estimates, (:rel_kt, :kep, :kep_rr))
    model(x, p) = model_referenceregion(conv; t=x, crr=crr,
        parameters=(rel_kt=p[1], kep=p[2], kep_rr=p[3]))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalues = [init_rel_kt[idx], init_kep[idx], init_kep_rr[idx]]
        rel_kt[idx], kep[idx], kep_rr[idx] = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return(estimates=(rel_kt=rel_kt, rel_ve=rel_ve, kep=kep, kep_rr=kep_rr), )
end

function fit_referenceregion_lls(; t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true)
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

function fit_constrained_referenceregion_nls(conv::Function=expconv;
    t::AbstractVector, crr::AbstractVector, ct::AbstractArray, mask=true, kep_rr=0.0)

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
        rrm_estimates = fit_referenceregion_nls(conv; t=t, crr=crr, ct=ct, mask=mask).estimates
        viable_estimates = positive_only_mask(rrm_estimates)
        kep_rr = interquartile_mean(rrm_estimates.kep_rr[viable_estimates])
    end
    lls_estimates = fit_constrained_referenceregion_lls(t=t, crr=crr, ct=ct, kep_rr=kep_rr, mask=mask).estimates
    init_rel_kt, init_kep = select(lls_estimates, (:rel_kt, :kep))
    model(x, p) = model_referenceregion(conv; t=x, crr=crr,
        parameters=(rel_kt=p[1], kep=p[2], kep_rr=kep_rr))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalues = [init_rel_kt[idx], init_kep[idx]]
        rel_kt[idx], kep[idx] = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return(estimates=(rel_kt=rel_kt, rel_ve=rel_ve, kep=kep, kep_rr=kep_rr), )
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
    return(estimates=(rel_kt=rel_kt, rel_ve=rel_ve, kep=kep, kep_rr=kep_rr), )
end


