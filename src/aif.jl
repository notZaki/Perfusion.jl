const _parker_parameters = (
    A = [0.809, 0.330],
    T = [0.17046, 0.365],
    σ = [0.0563, 0.132],
    α = 1.050,
    β = 0.1685,
    s = 38.078,
    τ = 0.483,
)

"""
    aif_parker(timepoints; parameters, hematocrit = 0.42)

Returns the concentration in blood plasma using the model defined in Parker et al. (2006, MRM 56(5). DOI: 10.1002/mrm.21066).
The `timepoints` input is in minutes with the bolus arriving at 0.
The model parameters from the paper will be used by default unless an optional `parameters` input is provided.
The `hematocrit` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_parker(timepoints; parameters = _parker_parameters, hematocrit = 0.42)
    @extract (A, T, σ, α, β, s, τ) parameters
    Cb = [(A[1] / (σ[1] * sqrt(2 * pi))) * exp(-(t - T[1])^2 / (2 * (σ[1])^2)) +
          (A[2] / (σ[2] * sqrt(2 * pi))) * exp(-(t - T[2])^2 / (2 * (σ[2])^2)) +
          α * exp(-β * t) / (1.0 + exp(-s * (t - τ))) for t in timepoints]
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

const _georgiou_parameters = (
    A = [0.37, 0.33, 10.06],
    m = [0.11, 1.17, 16.02],
    α = 5.26,
    β = 0.032,
    τ = 0.129,
)

"""
    aif_georgiou(timepoints; parameters, hematocrit = 0.35)

Returns the concentration in blood plasma using the model defined in Georgiou et al. (2019, MRM 81(3), DOI: 10.1002/mrm.27524).
The `timepoints` input is in minutes with the bolus arriving at 0.
The model parameters from the paper will be used by default unless an optional `parameters` input is provided.
The `hematocrit` is used to convert the concentration in blood to plasma and has default value of 0.35.
"""
function aif_georgiou(timepoints; parameters = _georgiou_parameters, hematocrit = 0.35)
    @extract (A, m, α, β, τ) parameters
    # The AIF can be split into two components or parts
    exponential_part = zeros(size(timepoints))
    gammavariate_part = zeros(size(timepoints))

    for (i, timepoint) in enumerate(timepoints)
        exponential_part[i] = sum(A .* exp.(-m * timepoint))
        if timepoint > 7.5
            # Gamma variate part undefined for high N (i.e. t > 7.5 min), so fix its value
            gammavariate_part[i] = 3.035
            continue
        end
        N = fld(timepoint, τ) # `fld` might be more efficient than `floor`
        for j = 0:N
            gammavariate_part[i] += gammavariate((j + 1) * α + j, β, timepoint - j * τ)
        end
    end
    Cb = exponential_part .* gammavariate_part
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

function gammavariate(alpha, beta, t)::Float64
    if t > 0
        gammavalue = (t^alpha * exp(-t / beta)) / (beta^(alpha + 1) * gamma(alpha + 1))
        if !isnan(gammavalue)
            return gammavalue
        end
    end
    return 0
end

"""
    aif_biexponential(timepoints; parameters, hematocrit = 0)

Returns an AIF defined by the biexponential model with the form ``D * (a_1 exp(-m_1 t) + a_2 exp(-m_2 t))``.
The `timepoints` input is in minutes with the bolus arriving at 0.
The `parameters` can be a named tuple, e.g. `parameters = (D=0.1, a=[4.0, 4.78], m=[0.144, 0.011])`.
The `hematocrit` is used to convert the concentration in blood to plasma.
The default value is set to to zero, i.e. no conversion takes place.
"""
function aif_biexponential(timepoints; parameters, hematocrit = 0)
    @extract (D, a, m) parameters
    Cb = [D * (a[1] * exp(-m[1] * t) + a[2] * exp(-m[2] * t)) for t in timepoints]
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

"""
    aif_weinmann(timepoints)

Returns an AIF defined by the biexponential model with parameters from
Tofts & Kermode (1991, MRM 17(2), DOI: 10.1002/mrm.1910170208) and
Weinmann et al. (1984, Phys. Chem. Phys. Med. NMR 16).
The `timepoints` input is in minutes with the bolus arriving at 0.
The model describes concentration in plasma, therefore a hematocrit is not necessary.
"""
function aif_weinmann(timepoints)
    weinmann_parameters = (; D = 0.1, a = [3.99, 4.78], m = [0.144, 0.0111])
    return aif_biexponential(timepoints; parameters = weinmann_parameters, hematocrit = 0)
end

"""
    aif_fritzhansen(timepoints)

Returns an AIF defined by the biexponential model with parameters based on data from
Fritz-Hansen et al. (1996, MRM 36(2), DOI: 10.1002/mrm.1910360209).
The model parameters were obtained from Whitcher & Schmid (2011, Journal of Statistical Software, 44(5))
The `timepoints` input is in minutes with the bolus arriving at 0.
The model describes concentration in plasma, therefore a hematocrit is not necessary.
"""
function aif_fritzhansen(timepoints; hematocrit = 0)
    fritzhansen_parameters = (; D = 1.0, a = [2.4, 0.62], m = [3.0, 0.016])
    return aif_biexponential(
        timepoints;
        parameters = fritzhansen_parameters,
        hematocrit = hematocrit,
    )
end

"""
    aif_orton1(timepoints; parameters, hematocrit = 0.42)

Returns the concentration in blood plasma using the bi-exponential model defined in
Orton et al. (2008, Phys. Med. Biol. 53(5), DOI: 10.1088/0031-9155/53/5/005).
The `timepoints` input is in minutes with the bolus arriving at 0.
The model parameters from the paper will be used by default unless an optional `parameters` input is provided.
The `hematocrit` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_orton1(
    timepoints;
    parameters = (aB = 10.4, μB = 12.0, aG = 1.24, μG = 0.169),
    hematocrit = 0.42,
)
    @extract (aB, μB, aG, μG) parameters
    AB = aB - aB * aG / (μB - μG)
    AG = aB * aG / (μB - μG)
    Cb = [AB * exp(-μB * t) + AG * exp(-μG * t) for t in timepoints]
    Cb[timepoints.<0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

"""
    aif_orton2(timepoints; parameters, hematocrit = 0.42)

Returns the concentration in blood plasma using the 2nd model defined in
Orton et al. (2008, Phys. Med. Biol. 53(5), DOI: 10.1088/0031-9155/53/5/005).
The `timepoints` input is in minutes with the bolus arriving at 0.
The model parameters from the paper will be used by default unless an optional `parameters` input is provided.
The `hematocrit` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_orton2(
    timepoints;
    parameters = (aB = 344, μB = 20.2, aG = 1.24, μG = 0.172),
    hematocrit = 0.42,
)
    @extract (aB, μB, aG, μG) parameters
    AB = aB - aB * aG / (μB - μG)
    AG = aB * aG / (μB - μG)^2
    Cb = [AB * t * exp(-μB * t) + AG * (exp(-μG * t) - exp(-μB * t)) for t in timepoints]
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

"""
    aif_orton3(timepoints; parameters, hematocrit = 0.42)

Returns the concentration in blood plasma using the 3rd model defined in
Orton et al. (2008, Phys. Med. Biol. 53(5), DOI: 10.1088/0031-9155/53/5/005).
The `timepoints` input is in minutes with the bolus arriving at 0.
The model parameters from the paper will be used by default unless an optional `parameters` input is provided.
The `hematocrit` is used to convert the concentration in blood to plasma and has default value of 0.42.
"""
function aif_orton3(
    timepoints;
    parameters = (aB = 2.84, μB = 22.8, aG = 1.36, μG = 0.171),
    hematocrit = 0.42,
)
    @extract (aB, μB, aG, μG) parameters
    tb = 2 * pi / μB
    Cb = zeros(length(timepoints))
    for (i, t) in enumerate(timepoints)
        if t <= tb
            Cb[i] = aB * (1 - cos(μB * t)) + aB * aG * _orton3_f(t, μG, μB)
        else
            Cb[i] = aB * aG * _orton3_f(tb, μG, μB) * exp(-μG * (t - tb))
        end
    end
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

function _orton3_f(t, α, μB)
    (1 / α) * (1 - exp(-α * t)) -
    (1 / (α^2 + μB^2)) * (α * cos(μB * t) + μB * sin(μB * t) - α * exp(-α * t))
end
