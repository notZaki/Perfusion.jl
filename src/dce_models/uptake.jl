function model_uptake(
    conv::Function = expconv;
    t::AbstractVector,
    parameters::NamedTuple,
    ca::AbstractVector,
)
    @extract (fp, ps, vp) parameters
    Tp = vp / (fp + ps)
    E = ps / (fp + ps)
    ca_conv_1 = zeros(length(t))
    for idx = 2:length(ca)
        ca_conv_1[idx] = ca_conv_1[idx-1] + (t[idx] - t[idx-1]) * ca[idx]
    end
    ct = fp .* ((1 - E) .* expconv(ca, 1 / Tp, t) .+ E .* ca_conv_1)
    return ct
end

function fit_uptake_nls(
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
    fp, ps, vp = (fill(NaN, volume_size...) for _ = 1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)
    model(x, p) =
        model_uptake(t = x, ca = ca, parameters = (fp = p[1], ps = p[2], vp = p[3]))
    lls_estimates = fit_uptake_lls(t = t, ca = ca, ct = ct, mask = mask).estimates
    init_fp, init_ps, init_vp = select(lls_estimates, (:fp, :ps, :vp))
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        initialvalue = [init_fp[idx], init_ps[idx], init_vp[idx]]
        (fp[idx], ps[idx], vp[idx]) = curve_fit(model, t, ct[idx, :], initialvalue).param
    end
    return (estimates = (fp = fp, ps = ps, vp = vp),)
end

function fit_uptake_lls(
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
    x1, fp, x3, ps, vp = (fill(NaN, volume_size...) for _ = 1:5)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 3)
    M[:, 2] = cumul_integrate(t, ca, TrapezoidalFast())
    M[:, 3] = cumul_integrate(t, M[:, 2], TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 1] = -cumul_integrate(t, ct[idx, :], TrapezoidalFast())
        (x1[idx], fp[idx], x3[idx]) = M \ ct[idx, :]
    end
    denum = @. x1 * fp - x3
    @. ps = fp * x3 / denum
    @. vp = fp^2 / denum
    @. ps[denum==0] = 0
    @. vp[denum==0] = 0
    return (estimates = (fp = fp, ps = ps, vp = vp),)
end
