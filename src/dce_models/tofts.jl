function model_tofts(
    conv::Function = expconv;
    t::AbstractVector,
    params::NamedTuple,
    ca::AbstractVector,
)
    if haskey(params, :ve)
        @extract (kt, ve) params
        kep = kt / ve
    else
        @extract (kt, kep) params
    end
    ct = _model_tofts(conv, t, [kt, kep], ca)
    if haskey(params, :vp)
        ct .+= params.vp * ca
    end
    return ct
end

function model_tofts2(
    t::AbstractVector,
    params::NamedTuple,
    ca::AbstractVector,
    conv::Function = expconv,
)
    if haskey(params, :ve)
        @extract (kt, ve) params
        kep = kt / ve
    else
        @extract (kt, kep) params
    end
    ct = _model_tofts(conv, t, [kt, kep], ca)
    if haskey(params, :vp)
        ct .+= params.vp * ca
    end
    return ct
end

function _model_tofts(
    conv::Function,
    t::AbstractVector,
    params::AbstractVector,
    ca::AbstractVector,
)
    kt, kep = params
    ct = kt * conv(ca, kep, t)
    return ct
end

function _model_extended_tofts(
    conv::Function,
    t::AbstractVector,
    params::AbstractVector,
    ca::AbstractVector,
)
    kt, kep, vp = params
    ct = kt * conv(ca, kep, t)
    ct .+= vp * ca
    return ct
end

function compute_ve_or_kep_if_it_is_missing!(collection)
    for (key, values) in collection
        if !haskey(values, :kep)
            collection[key] = (values..., kep = values.kt ./ values.ve)
        end
        if !haskey(values, :ve)
            collection[key] = (values..., ve = values.kt ./ values.kep)
        end
    end
    return collection
end

function compute_ve_or_kep_if_it_is_missing(namedtuple::NamedTuple)
    if !haskey(namedtuple, :kep)
        return (namedtuple..., kep = namedtuple.kt ./ namedtuple.ve)
    elseif !haskey(namedtuple, :ve)
        return (namedtuple..., ve = namedtuple.kt ./ namedtuple.kep)
    else
        return namedtuple
    end
end

function fit_tofts_nls(
    conv::Function = expconv;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) = resolve_fitting_inputs(; t, ca, ct, mask)
    kt, kep = (fill(NaN, volume_size...) for _ = 1:2)
    model(x, p) = _model_tofts(conv, x, p, ca)
    lls_est = fit_tofts_lls(; t, ca, ct, mask).est
    init_kt, init_kep = select(lls_est, (:kt, :kep))
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        initialvalues = [init_kt[idx], init_kep[idx]]
        (kt[idx], kep[idx]) = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    est = compute_ve_or_kep_if_it_is_missing((; kt, kep))
    return (; est)
end

function fit_extendedtofts_nls(
    conv::Function = expconv;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) = resolve_fitting_inputs(; t, ca, ct, mask)
    kt, kep, vp = (fill(NaN, volume_size...) for _ = 1:3)
    model(x, p) = _model_extended_tofts(conv, x, p, ca)
    lls_est = fit_extendedtofts_lls(; t, ca, ct, mask).est
    init_kt, init_kep, init_vp = select(lls_est, (:kt, :kep, :vp))
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        initialvalues = [init_kt[idx], init_kep[idx], init_vp[idx]]
        (kt[idx], kep[idx], vp[idx]) = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    est = compute_ve_or_kep_if_it_is_missing((; kt, kep, vp))
    return (; est)
end

function fit_tofts_lls(
    ;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) = resolve_fitting_inputs(; t, ca, ct, mask)
    kt, kep = (fill(NaN, volume_size...) for _ = 1:2)

    M = zeros(num_timepoints, 2)
    M[:, 1] = cumul_integrate(t, ca, TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        (kt[idx], kep[idx]) = M \ ct[idx, :]
    end
    est = compute_ve_or_kep_if_it_is_missing((; kt, kep))
    return (; est)
end

function fit_extendedtofts_lls(
    ;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) = resolve_fitting_inputs(; t, ca, ct, mask)
    kt, kep, vp = (fill(NaN, volume_size...) for _ = 1:3)

    M = zeros(num_timepoints, 3)
    M[:, 1] = cumul_integrate(t, ca, TrapezoidalFast())
    M[:, 3] = ca
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        (kt[idx], kep[idx], vp[idx]) = M \ ct[idx, :]
    end
    # Apply correction because fit actually returns: kt + kep * vp
    @. kt = kt - (kep * vp)
    est = compute_ve_or_kep_if_it_is_missing((; kt, kep, vp))
    return (; est)
end
