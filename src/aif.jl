const _parker_params = (
    A = [0.809, 0.330],
    T = [0.17046, 0.365],
    σ = [0.0563, 0.132],
    α = 1.050,
    β = 0.1685,
    s = 38.078,
    τ = 0.483,
)

"""
    aif_parker(t; params, hct = 0.42)

Returns the arterial input function defined in Parker et al. [@Parker2006].
The timepoints `t` are in minutes with the bolus arrival time defined as `t = 0`.
The model parameters from Ref [@Parker2006] will be used by default unless an optional `params` input is provided.
The hematocrit `hct` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_parker(t; params = _parker_params, hct = 0.42)
    @extract (A, T, σ, α, β, s, τ) params
    cb = [
        (A[1] / (σ[1] * sqrt(2 * pi))) * exp(-(t - T[1])^2 / (2 * (σ[1])^2)) +
        (A[2] / (σ[2] * sqrt(2 * pi))) * exp(-(t - T[2])^2 / (2 * (σ[2])^2)) +
        α * exp(-β * t) / (1.0 + exp(-s * (t - τ))) for t in t
    ]
    cb[t.<=0] .= 0
    ca = cb ./ (1 - hct)
    return ca
end

const _georgiou_params =
    (A = [0.37, 0.33, 10.06], m = [0.11, 1.17, 16.02], α = 5.26, β = 0.032, τ = 0.129)

"""
    aif_georgiou(t; params, hct = 0.35)

Returns the concentration in blood plasma using the model defined in Georgiou et al. [@Georgiou2018].
The timepoints `t` are in minutes with the bolus arrival defined as `t = 0`.
The model parameters from the paper will be used by default unless an optional `params` input is provided.
The hematocrit `hct` is used to convert the concentration in blood to plasma and has default value of 0.35.
"""
function aif_georgiou(t; params = _georgiou_params, hct = 0.35)
    @extract (A, m, α, β, τ) params
    # The AIF can be split into two parts
    exponential_part = zeros(size(t))
    gammavariate_part = zeros(size(t))

    for (i, t) in enumerate(t)
        exponential_part[i] = sum(A .* exp.(-m * t))
        if t > 7.5
            # Gamma variate part undefined for high N (i.e. t > 7.5 min), so fix its value
            gammavariate_part[i] = 3.035
            continue
        end
        N = fld(t, τ) # `fld` might be more efficient than `floor`
        for j = 0:N
            gammavariate_part[i] += gammavariate((j + 1) * α + j, β, t - j * τ)
        end
    end
    cb = exponential_part .* gammavariate_part
    cb[t.<=0] .= 0
    ca = cb ./ (1 - hct)
    return ca
end

function gammavariate(α, β, t)::Float64
    if t > 0
        gammavalue = (t^α * exp(-t / β)) / (β^(α + 1) * gamma(α + 1))
        if !isnan(gammavalue)
            return gammavalue
        end
    end
    return 0
end

"""
    aif_biexponential(t; params, hct = 0)

Returns an AIF defined by the biexponential model with the form ``D * (a_1 exp(-m_1 t) + a_2 exp(-m_2 t))``.
The `t` input is in minutes with the bolus arriving at 0.
The `params` can be a named tuple, e.g. `params = (D=0.1, a=[4.0, 4.78], m=[0.144, 0.011])`.
The `hct` is used to convert the concentration in blood to plasma. 
The default value is set to to zero, i.e. no conversion takes place.
"""
function aif_biexponential(t; params, hct = 0)
    @extract (D, a, m) params
    cb = [D * (a[1] * exp(-m[1] * t) + a[2] * exp(-m[2] * t)) for t in t]
    cb[t.<=0] .= 0
    ca = cb ./ (1 - hct)
    return ca
end

"""
    aif_weinmann(t)

Returns an AIF defined by the biexponential model with parameters from
Tofts & Kermode [@Tofts1991] and Weinmann et al. [@Weinmann1984].
The `t` input is in minutes with the bolus arriving at 0.
The model describes concentration in plasma, therefore a hematocrit is not necessary.
"""
function aif_weinmann(t)
    params = (; D = 0.1, a = [3.99, 4.78], m = [0.144, 0.0111])
    return aif_biexponential(t; params, hct = 0)
end

"""
    aif_fritzhansen(t)

Returns an AIF defined by the biexponential model with parameters based on data from
Fritz-Hansen et al. [@Fritz-Hansen1996].
The model parameters were adopted from Whitcher & Schmid [@Whitcher2011]
The `t` input is in minutes with the bolus arriving at 0.
The model describes concentration in plasma, therefore a hematocrit is not necessary.
"""
function aif_fritzhansen(t; hct = 0)
    params = (; D = 1.0, a = [2.4, 0.62], m = [3.0, 0.016])
    return aif_biexponential(t; params, hct)
end

"""
    aif_orton1(t; params = (aB = 10.4, μB = 12.0, aG = 1.24, μG = 0.169), hct = 0.42)

Returns the concentration in blood plasma using the bi-exponential model defined in
Orton et al. [@Orton2008].
The timepoints `t` are in minutes with the bolus arriving at 0.
The model parameters from the paper will be used by default unless an optional `params` input is provided.
The `hct` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_orton1(t; params = (aB = 10.4, μB = 12.0, aG = 1.24, μG = 0.169), hct = 0.42)
    @extract (aB, μB, aG, μG) params
    AB = aB - aB * aG / (μB - μG)
    AG = aB * aG / (μB - μG)
    cb = [AB * exp(-μB * t) + AG * exp(-μG * t) for t in t]
    cb[t.<0] .= 0
    ca = cb ./ (1 - hct)
    return ca
end

"""
    aif_orton2(t; params = (aB = 344, μB = 20.2, aG = 1.24, μG = 0.172), hct = 0.42)

Returns the concentration in blood plasma using the 2nd model defined in
Orton et al. [@Orton2008].
The timepoints `t` are in minutes with the bolus arriving at 0.
The model parameters from the paper will be used by default unless an optional `params` input is provided.
The `hematocrit` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_orton2(t; params = (aB = 344, μB = 20.2, aG = 1.24, μG = 0.172), hct = 0.42)
    @extract (aB, μB, aG, μG) params
    AB = aB - aB * aG / (μB - μG)
    AG = aB * aG / (μB - μG)^2
    cb = [AB * t * exp(-μB * t) + AG * (exp(-μG * t) - exp(-μB * t)) for t in t]
    cb[t.<=0] .= 0
    ca = cb ./ (1 - hct)
    return ca
end

"""
    aif_orton3(t; params = (aB = 2.84, μB = 22.8, aG = 1.36, μG = 0.171), hct = 0.42)

Returns the concentration in blood plasma using the 3rd model defined in
Orton et al. [@Orton2008].
The timepoints `t` are in minutes with the bolus arriving at `t = 0`.
The model parameters from the paper will be used by default unless an optional `params` input is provided.
The `hct` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_orton3(t; params = (aB = 2.84, μB = 22.8, aG = 1.36, μG = 0.171), hct = 0.42)
    @extract (aB, μB, aG, μG) params
    tb = 2 * pi / μB
    cb = zeros(length(t))
    for (i, t) in enumerate(t)
        if t <= tb
            cb[i] = aB * (1 - cos(μB * t)) + aB * aG * _orton3_f(t, μG, μB)
        else
            cb[i] = aB * aG * _orton3_f(tb, μG, μB) * exp(-μG * (t - tb))
        end
    end
    cb[t.<=0] .= 0
    ca = cb ./ (1 - hct)
    return ca
end

function _orton3_f(t, α, μB)
    (1 / α) * (1 - exp(-α * t)) -
    (1 / (α^2 + μB^2)) * (α * cos(μB * t) + μB * sin(μB * t) - α * exp(-α * t))
end
