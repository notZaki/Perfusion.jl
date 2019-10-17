function resolve_mask_size(mask, desired_size)
    if size(mask) == desired_size
        return mask .> 0
    elseif size(mask) == ()
        return repeat([mask .> 0], desired_size...)
    else
        error("Mask size: $(size(mask)) does not match input size $(desired_size)")
    end
end

function spgr(; M0, angle, TR, R1, R2star=0.0, TE=0.0)
    return @. M0 * sin(angle) * exp(-R2star*TE) * (1-exp(-R1*TR)) / (1-cos(angle)*exp(-R1*TR))
end

function conc_to_R1(Ct; r1, R10)
    return @. R10 + r1 * Ct
end

function R1_to_conc(R1; R10, r1)
    return @. (R1 - R10) / r1
end

function R1_to_signal(R1; S0, α, TR)
    return spgr(S0=S0, α=α, TR=TR, R1=R1)
end

function signal_to_R1(S; R10, α, TR, BAF::Int=1)
    (sX, sY, sZ, sT) = size(S)
    R1 = zeros(sX, sY, sZ, sT)
    S_pre = mean(S[:,:,:,1:BAF], dims=4)
    cosα = cos(α)
    for t=1:sT, k=1:sZ, j=1:sY, i=1:sX
        E1 = exp(-TR*R10[i,j,k])
        s = S[i,j,k,t] / S_pre[i,j,k]
        R1[i,j,k,t] = -(1/TR) * real(log( complex((1-s+s*E1-cosα*E1)/(1-s*cosα+s*E1*cosα-E1*cosα)) ))
    end
    return R1
end

function conc_to_signal(Ct; r1, S0, α, TR, R10)
    R1 = conc_to_R1(Ct; r1=r1, R10 = R10)
    signal = R1_to_signal(R1; S0=S0, α=α, TR=TR)
    return signal
end

function signal_to_conc(S; R10, α, TR, r1, BAF::Int=1)
    R1 = signal_to_R1(S; R10=R10, α=α, TR=TR, BAF=BAF)
    Ct = R1_to_conc(R1; R10=R10, r1=r1)
    return Ct
end

function fit_relaxation_despot(; signal::AbstractArray, angles::AbstractVector, TR::Real, mask=true)
    @assert length(angles) == size(signal)[end]
    num_angles = length(angles)
    if typeof(signal) <: AbstractVector
        @assert length(signal) == num_angles
        signal = reshape(signal, 1, num_angles)
    end
    volume_size = size(signal)[1:end-1]
    M0, T1 = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)

    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M0[idx], T1[idx] = _despot(signal[idx,:], angles, TR)
    end
    return(estimates=(M0=M0, T1=T1), dummy=0)
end

# Refs: Deoni, S. C. L., Peters, T. M., & Rutt, B. K. (2005). MRM 53(1), 237–241. https://doi.org/10.1002/mrm.20314
function _despot(signal::AbstractVector, alpha::AbstractVector, TR)
    x = signal ./ tan.(alpha)
    y = signal ./ sin.(alpha)
    numerator = sum(x .* y) - sum(x) .* sum(y) / length(alpha)
    denominator = sum(x.^2) - sum(x).^2 / length(alpha)  
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

function fit_relaxation_novifast(; signal::AbstractArray, angles::AbstractVector, TR::Real, mask=true)
    @assert length(angles) == size(signal)[end]
    num_angles = length(angles)
    if typeof(signal) <: AbstractVector
        @assert length(signal) == num_angles
        signal = reshape(signal, 1, num_angles)
    end
    volume_size = size(signal)[1:end-1]
    M0, T1 = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)

    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M0[idx], T1[idx] = _novifast(signal[idx,:], angles, TR)
    end
    return(estimates=(M0=M0, T1=T1), dummy=0)
end

# Refs: Ramos-Llorden, G., Vegas-Sanchez-Ferrero, G., Bjork, M., Vanhevel, F., Parizel, P. M., San Jose Estepar, R., … Sijbers, J. (2018). IEEE Transactions on Medical Imaging, 37(11), 2414–2427. https://doi.org/10.1109/TMI.2018.2833288
function _novifast(signal::AbstractVector, alpha::AbstractVector, TR; initialvalues=(M0=5000, T1=1500), maxiter=10, tol=1e-6)
    estimates = [
        initialvalues.M0 * (1 * exp(-TR / initialvalues.T1)), 
        exp(-TR/initialvalues.T1)
    ]
    sinα = sin.(alpha)
    cosα = cos.(alpha)
    k = 0
    relerr = 1e10
    while k <= maxiter && relerr > tol
        previous_estimates = copy(estimates)

        commondenom = 1 ./ (1 .- estimates[2] .* cosα);
        a = signal .* cosα .* commondenom;
        b = sinα .* commondenom;
        ahat = estimates[1] .* b .* cosα .* commondenom;
        svec = signal .* commondenom;

        svec_b = svec' * b;
        b_a = b' * a;
        svec_ahat = svec' * ahat;
        a_ahat = a' * ahat;
        b_b = b' * b;
        b_ahat = b' * ahat;

        num1 = svec_b*a_ahat - b_a*svec_ahat;
        num2 = b_b*svec_ahat - svec_b*b_ahat;
        denom = b_b*a_ahat - b_ahat*b_a;
        estimates = [num1/denom, num2/denom];

        k=k+1
        rel_err= norm(estimates .- previous_estimates) / norm(previous_estimates)
    end
    M0 = estimates[1] / (1 - estimates[2])
    if estimates[2] > 0
        T1 = -TR / log(estimates[2])
    else
        T1 = NaN
    end
    return M0, T1
end

function fit_relaxation_nls(;
    signal::AbstractArray,
    angles::AbstractVector,
    TR::Real,
    initialvalues=(M0=5000.0, T1=1500.0),
    mask=true
    )
    @assert length(angles) == size(signal)[end]
    num_angles = length(angles)
    if typeof(signal) <: AbstractVector
        @assert length(signal) == num_angles
        signal = reshape(signal, 1, num_angles)
    end
    volume_size = size(signal)[1:end-1]
    M0, T1 = (fill(NaN, volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)

    @. model(x,p) = spgr(M0=p[1], R1=1/p[2], angle=x, TR=TR)
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        fit = curve_fit(model, angles, signal[idx,:], [initialvalues.M0, initialvalues.T1])
        M0[idx], T1[idx] = fit.param 
    end
    return(estimates=(M0=M0, T1=T1), dummy=0)
end

const relaxation_dict = Dict{Symbol, Function}(
    :nls => fit_relaxation_nls,
    :despot => fit_relaxation_despot,
    :novifast => fit_relaxation_novifast
)

function fit_relaxation(fittype::Symbol; kwargs...)
    return relaxation_dict[fittype](;kwargs...)
end
