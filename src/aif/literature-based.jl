const _parker_parameters = (
    A = [0.809, 0.330],
    T = [0.17046, 0.365],
    σ = [0.0563, 0.132],
    α = 1.050,
    β = 0.1685,
    s = 38.078,
    τ = 0.483,
)

function aif_parker(timepoints; parameters = _parker_parameters, hematocrit = 0.42)
    @extract (A, T, σ, α, β, s, τ) parameters
    Cb = [( A[1] / (σ[1]*sqrt(2*pi)) ) * exp( -(t-T[1])^2 / (2*(σ[1])^2) ) +
          ( A[2] / (σ[2]*sqrt(2*pi)) ) * exp( -(t-T[2])^2 / (2*(σ[2])^2) ) +
          α * exp(-β*t) / (1. + exp( -s*(t-τ) )) for t in timepoints]
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

const _georgiou_parameters = (
    A = [0.37, 0.33, 10.06],
    m = [0.11, 1.17, 16.02],
    α = 5.26,
    β = 0.032,
    τ = 0.129
)

function aif_georgiou(timepoints; parameters = _georgiou_parameters, hematocrit = 0.35)
    @extract (A, m, α, β, τ) parameters
    #P = parameters; A = P.A; m = P.m; α = P.α; β = P.β; τ = P.τ;
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
        for j=0:N
            gammavariate_part[i] += gammavariate((j+1)*α+j, β, timepoint-j*τ)
        end
    end
    Cb = exponential_part .* gammavariate_part
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1-hematocrit)
    return Cp
end

function gammavariate(alpha, beta, t)::Float64
    if t>0
        gammavalue = (t^alpha * exp(-t/beta)) / (beta^(alpha + 1) * gamma(alpha + 1))
        if !isnan(gammavalue)
            return gammavalue
        end
    end
    return 0
end

function aif_biexponential(timepoints; parameters, hematocrit = 0.42)
    @extract (D, a, m) parameters
    Cb = [D * (a[1] * exp(-m[1]*t) + a[2] * exp(-m[2]*t)) for t in timepoints]
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

function aif_weinmann(timepoints; hematocrit = 0)
    weinmann_parameters = (
        D = 0.1,
        a = [3.99, 4.78],
        m = [0.144, 0.0111]
    )
    return aif_biexponential(timepoints, parameters = weinmann_parameters, hematocrit = hematocrit)
end

function aif_fritzhansen(timepoints; hematocrit = 0)
    fritzhansen_parameters = (
        D = 1.0,
        a = [2.4, 0.62],
        m = [3.0, 0.016]
    )
    return aif_biexponential(timepoints, parameters = fritzhansen_parameters, hematocrit = hematocrit)
end

function aif_orton1(timepoints;
    parameters = (aB = 10.4, μB = 12.0, aG =1.24, μG = 0.169),
    hematocrit = 0.42)
    @extract (aB, μB, aG, μG) parameters
    AB = aB - aB*aG / (μB-μG)
    AG = aB*aG / (μB-μG)
    Cb = [AB * exp(-μB*t) + AG * exp(-μG*t) for t in timepoints]
    Cb[timepoints.<0] .= 0
    Cp = Cb ./ ( 1 - hematocrit)
    return Cp
end

function aif_orton2(timepoints;
    parameters = (aB = 344, μB = 20.2, aG = 1.24, μG = 0.172),
    hematocrit = 0.42)
    @extract (aB, μB, aG, μG) parameters
    AB = aB - aB*aG / (μB-μG)
    AG = aB*aG / (μB-μG)^2
    Cb = [AB * t * exp(-μB*t) + AG * (exp(-μG*t) - exp(-μB*t)) for t in timepoints]
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

function aif_orton3(timepoints;
    parameters = (aB = 2.84, μB = 22.8, aG =1.36, μG = 0.171),
    hematocrit = 0.42)
    @extract (aB, μB, aG, μG) parameters
    tb = 2*pi / μB
    Cb = zeros(length(timepoints))
    for (i,t) in enumerate(timepoints)
        if t<=tb
            Cb[i] = aB * (1-cos(μB*t)) + aB*aG*_orton3_f(t,μG,μB)
        else
            Cb[i] = aB*aG*_orton3_f(tb, μG, μB)*exp(-μG*(t-tb))
        end
    end
    Cb[timepoints.<=0] .= 0
    Cp = Cb ./ (1 - hematocrit)
    return Cp
end

function _orton3_f(t, α, μB)
    (1/α) * (1-exp(-α*t)) - (1/(α^2 + μB^2)) * (α*cos(μB*t) + μB*sin(μB*t) - α*exp(-α*t))
end
