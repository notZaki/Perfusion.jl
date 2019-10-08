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
        f[i] = E*f[i-1] + A[i-1] * E0 + dA * E1
    end
    return f ./ B
end

function model_tofts(;t::AbstractVector, parameters::NamedTuple, Cp::AbstractVector)
    @extract (ktrans, kep) parameters
    Ct = ktrans * expconv(Cp, kep, t)
    if haskey(parameters, :vp)
        Ct .+= parameters.vp * Cp
    end
    return Ct
end

function fit_tofts(; t::AbstractVector, Cp::AbstractVector, Ct::AbstractArray, mask=true)
    @assert length(t) == length(Cp) == size(Ct)[end]
    num_timepoints = length(t)
    if typeof(Ct) <: AbstractVector
        @assert length(Ct) == num_timepoints
        Ct = reshape(Ct, 1, num_timepoints)
    end
    volume_size = size(Ct)[1:end-1]
    ktrans, kep = (zeros(volume_size...) for _=1:2)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 2)
    M[:,1] = cumul_integrate(t, Cp, TrapezoidalFast())
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, Ct[idx,:], TrapezoidalFast())
        (ktrans[idx], kep[idx]) = M \ Ct[idx,:]
    end
    return(estimates=(ktrans=ktrans, kep=kep), dummy=0)
end

function fit_extendedtofts(;t::AbstractVector, Cp::AbstractVector, Ct::AbstractArray, mask=true)
    @assert length(t) == length(Cp) == size(Ct)[end]
    num_timepoints = length(t)
    if typeof(Ct) <: AbstractVector
        @assert length(Ct) == num_timepoints
        Ct = reshape(Ct, 1, num_timepoints)
    end
    volume_size = size(Ct)[1:end-1]
    ktrans, kep, vp = (zeros(volume_size...) for _=1:3)
    resolved_mask = resolve_mask_size(mask, volume_size)

    M = zeros(num_timepoints, 3)
    M[:,1] = cumul_integrate(t, Cp, TrapezoidalFast())
    M[:,3] = Cp
    for idx in eachindex(IndexCartesian(), resolved_mask)
        if resolved_mask[idx] == false
            continue
        end
        M[:, 2] = -cumul_integrate(t, Ct[idx,:], TrapezoidalFast())
        (ktrans[idx], kep[idx], vp[idx]) = M \ Ct[idx,:]
    end
    # Apply correction because fit actually returns: ktrans + kep * vp
    @. ktrans = ktrans - (kep * vp)
    return(estimates=(ktrans=ktrans, kep=kep, vp=vp), dummy=0)
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
