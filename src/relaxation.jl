function spgr(; M0, angle, TR, R1, R2star = 0.0, TE = 0.0)
    return @. M0 * sin(angle) * exp(-R2star * TE) * (1 - exp(-R1 * TR)) /
              (1 - cos(angle) * exp(-R1 * TR))
end

function concentration_to_R1(Ct; r1, R10)
    return @. R10 + r1 * Ct
end

function R1_to_concentration(R1; R10, r1)
    return @. (R1 - R10) / r1
end

function R1_to_signal(R1; M0, angle, TR)
    return spgr(; M0, angle, TR, R1)
end

function signal_to_R1(signal; R10, angle, TR, BAF::Int = 1, mask = [true])
    if typeof(signal) <: AbstractVector
        signal = reshape(signal, 1, length(signal))
        input_was_vector = true
    else
        input_was_vector = false
    end
    volume_size = size(signal)[1:end-1]
    R1 = fill(NaN, size(signal)...)
    resolved_R10 = unify_size(R10, volume_size)
    resolved_mask = resolve_mask_size(mask, volume_size)

    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        R1[idx, :] = _signal_to_R1(signal[idx, :], resolved_R10[idx], angle, TR, BAF)
    end
    if input_was_vector
        R1 = R1[:]
    end
    return R1
end

function _signal_to_R1(signal, R10, angle, TR, BAF)
    S_pre = mean(signal[1:BAF])
    cosα = cos(angle)
    E1 = exp(-TR * R10)
    s = signal ./ S_pre
    log_term = @. (1 - s + s * E1 - cosα * E1) / (1 - s * cosα + s * E1 * cosα - E1 * cosα)
    R1 = [val > 0 ? (-1 / TR) * log(val) : NaN for val in log_term]
    return R1
end

function concentration_to_signal(Ct; r1, M0, angle, TR, R10, mask = true)
    R1 = concentration_to_R1(Ct; r1, R10, mask)
    signal = R1_to_signal(R1; M0, αngle, TR)
    return signal
end

function signal_to_concentration(S; R10, angle, TR, r1, BAF::Int = 1, mask = true)
    R1 = signal_to_R1(S; R10, angle, TR, BAF, mask)
    Ct = R1_to_concentration(R1; R10, r1)
    return Ct
end

function fit_relaxation_despot(
    ;
    signal::AbstractArray,
    angles::AbstractVector,
    TR::Real,
    mask = true,
)
    (signal, angles, mask, num_angles, volume_size) = resolve_relaxation_inputs(; signal, angles, mask)
    M0, T1 = (fill(NaN, volume_size...) for _ = 1:2)

    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        M0[idx], T1[idx] = _despot(signal[idx, :], angles, TR)
    end
    return (; estimates = (; M0, T1))
end

# Refs: Deoni, S. C. L., Peters, T. M., & Rutt, B. K. (2005). MRM 53(1), 237–241. https://doi.org/10.1002/mrm.20314
function _despot(signal::AbstractVector, alpha::AbstractVector, TR)
    x = signal ./ tan.(alpha)
    y = signal ./ sin.(alpha)
    numerator = sum(x .* y) - sum(x) .* sum(y) / length(alpha)
    denominator = sum(x .^ 2) - sum(x) .^ 2 / length(alpha)
    slope = numerator ./ denominator
    intercept = mean(y) - slope * mean(x)
    M0 = intercept / (1 - slope)
    if slope > 0
        T1 = -TR / log(slope)
    else
        T1 = NaN
    end
    return M0, T1
end

function fit_relaxation_novifast(
    ;
    signal::AbstractArray,
    angles::AbstractVector,
    TR::Real,
    mask = true,
)
    (signal, angles, mask, num_angles, volume_size) = resolve_relaxation_inputs(; signal, angles, mask)
    M0, T1 = (fill(NaN, volume_size...) for _ = 1:2)

    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        M0[idx], T1[idx] = _novifast(signal[idx, :], angles, TR)
    end
    return (; estimates = (; M0, T1))
end

# Refs: Ramos-Llorden, G., Vegas-Sanchez-Ferrero, G., Bjork, M., Vanhevel, F., Parizel, P. M., San Jose Estepar, R., … Sijbers, J. (2018). IEEE Transactions on Medical Imaging, 37(11), 2414–2427. https://doi.org/10.1109/TMI.2018.2833288
function _novifast(
    signal::AbstractVector,
    alpha::AbstractVector,
    TR;
    initialvalues = (M0 = 5000, T1 = 1500),
    maxiter = 10,
    tol = 1e-6,
)
    estimates = [
        initialvalues.M0 * (1 * exp(-TR / initialvalues.T1)),
        exp(-TR / initialvalues.T1),
    ]
    sinα = sin.(alpha)
    cosα = cos.(alpha)
    k = 0
    relerr = 1e10
    while k <= maxiter && relerr > tol
        previous_estimates = copy(estimates)

        commondenom = 1 ./ (1 .- estimates[2] .* cosα)
        a = signal .* cosα .* commondenom
        b = sinα .* commondenom
        ahat = estimates[1] .* b .* cosα .* commondenom
        svec = signal .* commondenom

        svec_b = svec' * b
        b_a = b' * a
        svec_ahat = svec' * ahat
        a_ahat = a' * ahat
        b_b = b' * b
        b_ahat = b' * ahat

        num1 = svec_b * a_ahat - b_a * svec_ahat
        num2 = b_b * svec_ahat - svec_b * b_ahat
        denom = b_b * a_ahat - b_ahat * b_a
        estimates = [num1 / denom, num2 / denom]

        k = k + 1
        rel_err = norm(estimates .- previous_estimates) / norm(previous_estimates)
    end
    M0 = estimates[1] / (1 - estimates[2])
    if estimates[2] > 0
        T1 = -TR / log(estimates[2])
    else
        T1 = NaN
    end
    return M0, T1
end


function fit_relaxation_nls(
    ;
    signal::AbstractArray,
    angles::AbstractVector,
    TR::Real,
    initialvalues = (M0 = 5000.0, T1 = 1500.0),
    mask = true,
)
    (signal, angles, mask, num_angles, volume_size) = resolve_relaxation_inputs(; signal, angles, mask)
    M0, T1 = (fill(NaN, volume_size...) for _ = 1:2)

    @. model(x, p) = _spgr(x, p[1], 1 / p[2], TR)
    for idx in eachindex(IndexCartesian(), mask)
        if mask[idx] == false
            continue
        end
        fit = curve_fit(model, angles, signal[idx, :], [initialvalues.M0, initialvalues.T1])
        M0[idx], T1[idx] = fit.param
    end
    return (; estimates = (; M0, T1))
end

# [Todo] Replace with spgr() --- but check for speed difference first
function _spgr(angle, M0, R1, TR)
    return @. M0 * sin(angle) * (1 - exp(-R1 * TR)) / (1 - cos(angle) * exp(-R1 * TR))
end

const relaxation_dict = Dict{Symbol,Function}(
    :nls => fit_relaxation_nls,
    :despot => fit_relaxation_despot,
    :novifast => fit_relaxation_novifast,
)

function fit_relaxation(fittype::Symbol = :novifast; kwargs...)
    return relaxation_dict[fittype](; kwargs...)
end
