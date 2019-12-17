function model_exchange(
    conv::Function = expconv;
    t::AbstractVector,
    parameters::NamedTuple,
    ca::AbstractVector,
)
    if all(haskey(parameters, key) for key in (:Tp, :Te))
        @extract (ve, vp) parameters
        fp = vp / parameters.Tp
        ps = ve / parameters.Te
    else
        @extract (fp, ps, ve, vp) parameters
    end
    ct = _model_exchange(conv, t, [fp, ps, ve, vp], ca)
    return ct
end

function _model_exchange(conv, t, parameters, ca)
    fp, ps, ve, vp = parameters
    Tp = vp / fp
    Te = ve / ps
    T = (vp + ve) / fp
    Tplus = (T + Te + sqrt((T + Te)^2 - 4 * Tp * Te)) / 2
    Tminus = (T + Te - sqrt((T + Te)^2 - 4 * Tp * Te)) / 2
    ct = ((T - Tminus) / (Tplus - Tminus)) .* conv(ca, 1 / Tplus, t)
    ct .+= ((Tplus - T) / (Tplus - Tminus)) .* conv(ca, 1 / Tminus, t)
    ct .*= fp
    return ct
end

function fit_exchange_nls(
    conv::Function = expconv;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    fp, ps, vp, ve = (fill(NaN, volume_size...) for _ = 1:4)
    model(x, p) = _model_exchange(conv, x, p, ca)
    lls_estimates = fit_exchange_lls(t = t, ca = ca, ct = ct, mask = mask).estimates
    init_fp, init_ps, init_ve, init_vp = select(lls_estimates, (:fp, :ps, :ve, :vp))
    resolved_mask = resolve_mask_size(mask, volume_size)
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalue = [init_fp[idx], init_ps[idx], init_ve[idx], init_vp[idx]]
        fp[idx], ps[idx], ve[idx], vp[idx] = curve_fit(model, t, ct[idx, :], initialvalue).param
    end
    T = @. (vp + ve) / fp
    Tp = @. vp / fp
    Te = @. ve / ps
    return (estimates = (fp = fp, ps = ps, ve = ve, vp = vp, T = T, Tp = Tp, Te = Te),)
end

function fit_exchange_lls(
    ;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    @assert length(t) == length(ca) == size(ct)[end]
    num_timepoints = length(t)
    if typeof(ct) <: AbstractVector
        @assert length(ct) == num_timepoints
        ct = reshape(ct, 1, num_timepoints)
    end
    volume_size = size(ct)[1:end-1]
    α, β, γ, fp = (fill(NaN, volume_size...) for _ = 1:4)
    resolved_mask = resolve_mask_size(mask, volume_size)
    M = zeros(num_timepoints, 4)
    M[:, 4] .= cumul_integrate(t, ca, TrapezoidalFast())
    M[:, 3] .= cumul_integrate(t, M[:, 4], TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 2] .= -cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        M[:, 1] .= cumul_integrate(t, M[:, 2], TrapezoidalFast())
        α[idx], β[idx], γ[idx], fp[idx] = M \ ct[idx, :]
    end
    T = @. γ / (α * fp)
    Te = @. (β / α) - T
    Tp = @. 1 / (α * Te)
    vp = @. fp * Tp
    ve = @. fp * (T - Tp)
    ps = @. ve / Te
    return (estimates = (fp = fp, ps = ps, ve = ve, vp = vp, T = T, Te = Te, Tp = Tp),)
end
