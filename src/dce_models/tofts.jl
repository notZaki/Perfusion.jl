function model_tofts(conv::Function=expconv; t::AbstractVector, parameters::NamedTuple, cp::AbstractVector)
    if haskey(parameters, :ve)
        @extract (kt, ve) parameters
        kep = kt / ve
    else
        @extract (kt, kep) parameters
    end
    ct = _model_tofts(conv, t, [kt, kep], cp)
    if haskey(parameters, :vp)
        ct .+= parameters.vp * cp
    end
    return ct
end

function _model_tofts(conv::Function, t::AbstractVector, parameters::AbstractVector, cp::AbstractVector)
    kt, kep = parameters
    ct = kt * conv(cp, kep, t)
    return ct
end

function _model_extended_tofts(conv::Function, t::AbstractVector, parameters::AbstractVector, cp::AbstractVector)
    kt, kep, vp = parameters
    ct = kt * conv(cp, kep, t)
    ct .+= vp * cp
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

function fit_tofts_nls(conv::Function=expconv;
    t::AbstractVector, cp::AbstractVector, ct::AbstractArray, mask=true)

    @assert length(t) == length(cp) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    kt, kep = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)
    model(x, p) = _model_tofts(conv, x, p, cp)
    lls_estimates = fit_tofts_lls(t=t, cp=cp, ct=ct, mask=mask).estimates
    init_kt, init_kep = select(lls_estimates, (:kt, :kep))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalues = [init_kt[idx], init_kep[idx]]
        (kt[idx], kep[idx]) = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    estimates = compute_ve_or_kep_if_it_is_missing((kt=kt, kep=kep))
    return(estimates=estimates, )
end

function fit_extendedtofts_nls(conv::Function=expconv;
    t::AbstractVector, cp::AbstractVector, ct::AbstractArray, mask=true)

    @assert length(t) == length(cp) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    kt, kep, vp = (fill(NaN, volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)

    model(x, p) = _model_extended_tofts(conv, x, p, cp)
    lls_estimates = fit_extendedtofts_lls(t=t, cp=cp, ct=ct, mask=mask).estimates
    init_kt, init_kep, init_vp = select(lls_estimates, (:kt, :kep, :vp))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalues = [init_kt[idx], init_kep[idx], init_vp[idx]]
        (kt[idx], kep[idx], vp[idx]) = curve_fit(model, t, ct[idx, :], initialvalues).param
    end
    estimates = compute_ve_or_kep_if_it_is_missing((kt=kt, kep=kep, vp=vp))
    return(estimates=estimates, )
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
    estimates = compute_ve_or_kep_if_it_is_missing((kt=kt, kep=kep))
    return(estimates=estimates, )
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
    estimates = compute_ve_or_kep_if_it_is_missing((kt=kt, kep=kep, vp=vp))
    return(estimates=estimates, )
end

