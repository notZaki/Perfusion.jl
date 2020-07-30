function model_filtration(
    conv::Function = expconv;
    t::AbstractVector,
    params::NamedTuple,
    ca::AbstractVector,
)
    if all(haskey(params, key) for key in (:Tp, :Te))
        @extract (ve, vp) params
        fp = vp / params.Tp
        ps = ve / params.Te
    else
        @extract (fp, ps, ve, vp) params
    end
    ct = _model_filtration(conv, t, [fp, ps, ve, vp], ca)
    return ct
end

function _model_filtration(conv, t, params, ca)
    fp, ps, ve, vp = params
    Tminus = vp / fp
    Tplus = ve / ps
    T = (vp + ve) / fp
    ct = ((T - Tminus) / (Tplus - Tminus)) .* expconv(ca, 1 / Tplus, t)
    ct .+= ((Tplus - T) / (Tplus - Tminus)) .* expconv(ca, 1 / Tminus, t)
    ct .*= fp
    return ct
end

function fit_filtration_nls(
    conv::Function = expconv;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) = resolve_fitting_inputs(; t, ca, ct, mask)
    fp, ps, vp, ve = (fill(NaN, volume_size...) for _ = 1:4)
    model(x, p) = _model_filtration(conv, x, p, ca)
    lls_est = fit_filtration_lls(; t, ca, ct, mask).est
    init_fp, init_ps, init_ve, init_vp = select(lls_est, (:fp, :ps, :ve, :vp))
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        initialvalue = [init_fp[idx], init_ps[idx], init_ve[idx], init_vp[idx]]
        fp[idx], ps[idx], ve[idx], vp[idx] = curve_fit(model, t, ct[idx, :], initialvalue).param
    end
    T = @. (vp + ve) / fp
    Tp = @. vp / fp
    Te = @. ve / ps
    @. T[fp==0] = 0
    @. Tp[fp==0] = 0
    @. Te[ps==0] = 0
    return (; est = (; fp, ps, ve, vp, T, Tp, Te))
end

function fit_filtration_lls(
    ;
    t::AbstractVector,
    ca::AbstractVector,
    ct::AbstractArray,
    mask = true,
)
    (t, ct, mask, num_timepoints, volume_size) = resolve_fitting_inputs(; t, ca, ct, mask)
    α, β, γ, fp = (fill(NaN, volume_size...) for _ = 1:4)

    M = zeros(num_timepoints, 4)
    M[:, 4] .= cumul_integrate(t, ca, TrapezoidalFast())
    M[:, 3] .= cumul_integrate(t, M[:, 4], TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        M[:, 2] .= -cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        M[:, 1] .= cumul_integrate(t, M[:, 2], TrapezoidalFast())
        α[idx], β[idx], γ[idx], fp[idx] = M \ ct[idx, :]
    end
    root_term = @. β^2 - 4 * α
    @. root_term[root_term<0] = 0
    Tp = @. (β - sqrt(root_term)) / (2 * α)
    Te = @. (β + sqrt(root_term)) / (2 * α)
    T = @. γ / (α * fp)
    vp = @. fp * Tp
    ve = @. fp * (T - Tp)
    ps = @. ve / Te
    return (; est = (; fp, ps, ve, vp, T, Te, Tp))
end
