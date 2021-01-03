function model_referenceregion(
    conv::Function = expconv;
    t::AbstractVector,
    params::NamedTuple,
    crr::AbstractVector,
)
    @extract (rel_kt, kep, kep_rr) params
    ct = rel_kt .* crr .+ (rel_kt * (kep_rr - kep)) .* conv(crr, kep, t)
    return ct
end

function fit_rrm_nls(
    conv::Function = expconv;
    t::AbstractVector,
    crr::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) =
        resolve_fitting_inputs(; t, ca = crr, ct, mask)
    rel_kt, kep, kep_rr = (fill(NaN, volume_size...) for _ = 1:3)
    lls_est = fit_rrm_lls(; t, crr, ct, mask).est
    init_rel_kt, init_kep, init_kep_rr = select(lls_est, (:rel_kt, :kep, :kep_rr))
    model(x, p) = model_referenceregion(
        conv;
        t = x,
        crr,
        params = (rel_kt = p[1], kep = p[2], kep_rr = p[3]),
    )
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        initialvalues = [init_rel_kt[idx], init_kep[idx], init_kep_rr[idx]]
        rel_kt[idx], kep[idx], kep_rr[idx] =
            curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return (; est = (; rel_kt, rel_ve, kep, kep_rr))
end

function fit_rrm_lls(;
    t::AbstractVector,
    crr::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) =
        resolve_fitting_inputs(; t, ca = crr, ct, mask)
    x1, x2, x3 = (fill(NaN, volume_size...) for _ = 1:3)

    M = zeros(num_timepoints, 3)
    M[:, 1] = crr
    M[:, 2] = cumul_integrate(t, crr, TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        M[:, 3] = -cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        x1[idx], x2[idx], x3[idx] = M \ ct[idx, :]
    end
    # Apply correction because fit actually returns: [rel_kt, kt/ve_rr, kep]
    rel_kt = x1
    rel_ve = x2 ./ x3
    kep = x3
    kep_rr = x2 ./ x1
    return (; est = (; rel_kt, rel_ve, kep, kep_rr))
end

function fit_crrm_nls(
    conv::Function = expconv;
    t::AbstractVector,
    crr::AbstractVector,
    ct::AbstractArray,
    mask = true,
    kep_rr = 0.0,
)
    (t, ct, mask, num_timepoints, volume_size) =
        resolve_fitting_inputs(; t, ca = crr, ct, mask)
    rel_kt, kep = (fill(NaN, volume_size...) for _ = 1:2)

    if kep_rr <= 0
        rrm_est = fit_rrm_nls(conv; t, crr, ct, mask).est
        viable_est = positive_only_mask(rrm_est)
        kep_rr = interquartile_mean(rrm_est.kep_rr[viable_est])
    end
    lls_est = fit_crrm_lls(; t, crr, ct, kep_rr, mask = mask).est
    init_rel_kt, init_kep = select(lls_est, (:rel_kt, :kep))
    model(x, p) = model_referenceregion(
        conv;
        t = x,
        crr,
        params = (rel_kt = p[1], kep = p[2], kep_rr = kep_rr),
    )
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        initialvalues = [init_rel_kt[idx], init_kep[idx]]
        rel_kt[idx], kep[idx] = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return (; est = (; rel_kt, rel_ve, kep, kep_rr))
end

function fit_crrm_lls(;
    t::AbstractVector,
    crr::AbstractVector,
    ct::AbstractArray,
    mask = true,
    kep_rr = 0.0,
)
    (t, ct, mask, num_timepoints, volume_size) =
        resolve_fitting_inputs(; t, ca = crr, ct, mask)
    rel_kt, kep = (fill(NaN, volume_size...) for _ = 1:2)

    if kep_rr <= 0
        rrm_est = fit_rrm_lls(; t, crr, ct, mask).est
        viable_est = positive_only_mask(rrm_est)
        kep_rr = interquartile_mean(rrm_est.kep_rr[viable_est])
    end

    M = zeros(num_timepoints, 2)
    M[:, 1] = crr + kep_rr * cumul_integrate(t, crr, TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        rel_kt[idx], kep[idx] = M \ ct[idx, :]
    end
    rel_ve = @. rel_kt * kep_rr / kep
    return (est = (rel_kt = rel_kt, rel_ve = rel_ve, kep = kep, kep_rr = kep_rr),)
end

function fit_errm_lls(;
    t::AbstractVector,
    crr::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) =
        resolve_fitting_inputs(; t, ca = crr, ct, mask)
    x1, x2, x3, x4 = (fill(NaN, volume_size...) for _ = 1:4)

    M = zeros(num_timepoints, 4)
    M[:, 1] = cumul_integrate(t, crr, TrapezoidalFast())
    M[:, 2] = cumul_integrate(t, M[:, 1], TrapezoidalFast())
    M[:, 4] = crr
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        cur_ct = cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        M[:, 3] = -cumul_integrate(t, cur_ct, TrapezoidalFast())
        x1[idx], x2[idx], x3[idx], x4[idx] = M \ cur_ct
    end
    # Apply correction because fit actually returns: [mess, mess, kep, vp/kt_rr]
    a = x1 ./ x4
    b = x2 ./ x4
    root_term = @. a^2 - 4 * b
    @. root_term[root_term<0] = 0
    kep_rr = @. (a - sqrt(root_term)) / 2
    rel_kt = @. x4 * (a - x3 - kep_rr)
    rel_ve = @. rel_kt * kep_rr / x3
    kep = x3
    rel_vp = x4
    return (; est = (; rel_kt, rel_ve, rel_vp, kep, kep_rr))
end

function fit_cerrm_lls(;
    t::AbstractVector,
    crr::AbstractVector,
    ct::AbstractArray,
    mask = true,
    kep_rr = 0.0,
)
    (t, ct, mask, num_timepoints, volume_size) =
        resolve_fitting_inputs(; t, ca = crr, ct, mask)
    x1, rel_vp, kep = (fill(NaN, volume_size...) for _ = 1:3)

    if kep_rr <= 0
        rrm_est = fit_errm_lls(; t, crr, ct, mask).est
        viable_est = positive_only_mask(rrm_est)
        kep_rr = interquartile_mean(rrm_est.kep_rr[viable_est])
    end

    M = zeros(num_timepoints, 3)
    crr_int = cumul_integrate(t, crr, TrapezoidalFast())
    crr_int2 = cumul_integrate(t, crr_int, TrapezoidalFast())
    M[:, 1] = crr_int + kep_rr * crr_int2
    M[:, 2] = crr + kep_rr * crr_int
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        cur_ct = cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        M[:, 3] = -cumul_integrate(t, cur_ct, TrapezoidalFast())
        x1[idx], rel_vp[idx], kep[idx] = M \ cur_ct
    end
    rel_kt = @. x1 - kep * rel_vp
    rel_ve = @. rel_kt * kep_rr / kep
    return (; est = (; rel_kt, rel_ve, rel_vp, kep, kep_rr))
end

function fit_rrift_with_cerrm(; crr, cp, t, ct, tail_start, kep_rr = 0.0, mask = true)
    cerrm = fit_cerrm_lls(; crr, ct, t, kep_rr, mask).est
    kep_rr = cerrm.kep_rr
    kt_rr = fit_rrift(; t, cp, crr, kep_rr, tail_start)
    ve_rr = kt_rr / kep_rr
    est = relative_to_absolute(cerrm; kt_rr, ve_rr)
    return (est = (; est..., kep_rr, kt_rr, ve_rr))
end

function fit_rrift_with_crrm(; crr, cp, t, ct, tail_start, kep_rr = 0.0, mask = true)
    crrm = fit_crrm_lls(; crr, ct, t, kep_rr, mask).est
    kep_rr = crrm.kep_rr
    kt_rr = fit_rrift(; t, cp, crr, kep_rr, tail_start)
    ve_rr = kt_rr / kep_rr
    est = relative_to_absolute(crrm; kt_rr, ve_rr)
    return (est = (; est..., kep_rr, kt_rr, ve_rr))
end

function fit_rrift(;
    crr::AbstractVector,
    cp::AbstractVector,
    t::AbstractVector,
    tail_start::Int,
    kep_rr::Number,
)
    @assert length(crr) == length(cp) == length(t)
    @assert tail_start < length(crr)
    crr_tail = crr[tail_start:end]
    cp_tail = cp[tail_start:end]
    t_tail = t[tail_start:end]
    numerator =
        crr_tail .- crr_tail[1] .+
        kep_rr .* cumul_integrate(t_tail, crr_tail, TrapezoidalFast())
    denominator = cumul_integrate(t_tail, cp_tail, TrapezoidalFast())
    kt_rr = denominator \ numerator
    return kt_rr
end

function relative_to_absolute(rel_params; kt_rr, ve_rr)
    kt = @. rel_params.rel_kt * kt_rr
    ve = @. rel_params.rel_ve * ve_rr
    kep = rel_params.kep
    if haskey(rel_params, :rel_vp)
        vp = @. rel_params.rel_vp * kt_rr
    else
        vp = zeros(size(kt))
    end
    return (kt = kt, ve = ve, vp = vp, kep = kep)
end
