using Perfusion
using Statistics
using Test
using UnicodePlots

function is_all_nan(x::NamedTuple)
    for val in x
        if any(@. !isnan(val))
            return false
        end
    end
    return true
end

@testset "AIFs" begin
    scan_timepoints = collect(1:600) ./ 60
    aif_timepoints = scan_timepoints .- 1.0 # Define bolus arrival at 1 minute

    aifs_to_test = (
        "Fritz-Hansen" => aif_fritzhansen,
        "Georgiou" => aif_georgiou,
        "Orton V1" => aif_orton1,
        "Orton V2" => aif_orton2,
        "Orton V3" => aif_orton3,
        "Parker" => aif_parker,
        "Weinmann" => aif_weinmann,
    )
    plot_cfg(name) =
        (title = name, xlabel = "Gd [mM]", ylabel = "time [min]", canvas = DotCanvas)
    printdiv() = printstyled(
        "\n==================================================================================\n",
        color = :light_black,
    )

    for (name, func) in aifs_to_test
        println(lineplot(scan_timepoints, func(aif_timepoints); plot_cfg(name)...))
        printdiv()
    end
end

@testset "T1 fitting" begin
    T1 = 2500
    R1 = 1 / T1
    M0 = 10_000
    α = [5, 10, 15, 20, 25, 30]
    TR = 5

    signal = Perfusion.spgr(M0 = M0, R1 = R1, angle = α, TR = TR)

    despot = fit_relaxation(:despot, signal = signal, angles = α, TR = TR)
    novifast = fit_relaxation(:novifast, signal = signal, angles = α, TR = TR)
    nls = fit_relaxation(:nls, signal = signal, angles = α, TR = TR)

    @test round(despot.est.M0[1], digits = 3) == M0
    @test round(despot.est.T1[1], digits = 3) == T1

    @test round(novifast.est.M0[1], digits = 3) == M0
    @test round(novifast.est.T1[1], digits = 3) == T1

    @test round(nls.est.M0[1], digits = 3) == M0
    @test round(nls.est.T1[1], digits = 3) == T1
end

@testset "Signal <-> Concentration" begin
    scan_timepoints = collect(1:600) ./ 60
    aif_timepoints = scan_timepoints .- 1.0 # Define bolus arrival at 1 minute
    r1 = 3.3 / 1000 # units: mM/ms
    R10 = 1 / 1000 # units: 1/ms
    TR = 5 # units: ms
    angle = deg2rad(20)
    M0 = 1_000
    C = aif_georgiou(aif_timepoints)
    signal = concentration_to_signal(C; r1, angle, TR, R10, M0)
    new_C = signal_to_concentration(signal; r1, angle, TR, R10)
    @test isapprox(C, new_C)
end

@testset "Tofts model" begin
    t = collect(1:600) ./ 60
    ca = aif_parker(t .- 1)

    params = (kt = 0.25, kep = 0.5)
    ct = model_tofts(; t, ca, params)

    est = fit_model(:tofts, :lls; t, ca, ct, mask = [true]).est
    @test round(est.kt[1], digits = 3) == params.kt
    @test round(est.kep[1], digits = 3) == params.kep

    est = fit_model(:tofts; t, ca, ct).est
    @test round(est.kt[1], digits = 3) == params.kt
    @test round(est.kep[1], digits = 3) == params.kep

    @test is_all_nan(fit_model(:tofts, :lls; t, ca, ct, mask = false).est)
    @test is_all_nan(fit_model(:tofts, :nls; t, ca, ct, mask = [false]).est)
    @test_throws ErrorException fit_model(:tofts; t, ca, ct, mask = [true, true])
end

@testset "Extended Tofts model" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    params = (kt = 0.5, kep = 1.0, vp = 0.1)
    ct = model_tofts(; t, ca, params)

    est = fit_model(:extendedtofts, :lls; t, ca, ct, mask = [true]).est
    @test round(est.kt[1], digits = 3) == params.kt
    @test round(est.kep[1], digits = 3) == params.kep
    @test round(est.vp[1], digits = 3) == params.vp

    est = fit_model(:extendedtofts; t, ca, ct).est
    @test round(est.kt[1], digits = 3) == params.kt
    @test round(est.kep[1], digits = 3) == params.kep
    @test round(est.vp[1], digits = 3) == params.vp

    @test is_all_nan(fit_model(:extendedtofts, :lls; t, ca, ct, mask = false).est)
    @test is_all_nan(fit_model(:extendedtofts, :nls; t, ca, ct, mask = [false]).est)
    @test_throws ErrorException fit_model(:extendedtofts; t, ca, ct, mask = [true, true])
end

@testset "Compartmental tissue uptake model" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    params = (fp = 0.75, ps = 0.05, vp = 0.25)
    ct = model_uptake(; t, ca, params)

    est = fit_model(:uptake, :lls; t, ca, ct, mask = [true]).est
    @test round(est.fp[1], digits = 2) == params.fp
    @test round(est.ps[1], digits = 2) == params.ps
    @test round(est.vp[1], digits = 2) == params.vp

    est = fit_model(:uptake, :nls; t, ca, ct, mask = [true]).est
    @test round(est.fp[1], digits = 3) == params.fp
    @test round(est.ps[1], digits = 3) == params.ps
    @test round(est.vp[1], digits = 3) == params.vp

    @test is_all_nan(fit_model(:uptake, :lls; t, ca, ct, mask = false).est)
    @test is_all_nan(fit_model(:uptake, :nls; t, ca, ct, mask = [false]).est)
    @test_throws ErrorException fit_model(:uptake; t, ca, ct, mask = [true, true])
end

@testset "Two compartment exchange model" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    params = (fp = 0.75, ps = 0.05, vp = 0.25, ve = 0.10)
    ct = model_exchange(; t, ca, params)

    est = fit_model(:exchange, :lls; t, ca, ct, mask = [true]).est
    @test round(est.fp[1], digits = 3) == params.fp
    @test round(est.ps[1], digits = 3) == params.ps
    @test round(est.ve[1], digits = 3) == params.ve
    @test round(est.vp[1], digits = 3) == params.vp

    est = fit_model(:exchange, :nls; t, ca, ct, mask = [true]).est
    @test round(est.fp[1], digits = 3) == params.fp
    @test round(est.ps[1], digits = 3) == params.ps
    @test round(est.ve[1], digits = 3) == params.ve
    @test round(est.vp[1], digits = 3) == params.vp

    @test is_all_nan(fit_model(:exchange, :lls; t, ca, ct, mask = false).est)
    @test is_all_nan(fit_model(:exchange, :nls; t, ca, ct, mask = [false]).est)
    @test_throws ErrorException fit_model(:exchange; t, ca, ct, mask = [true, true])

    fp, ps, vp, ve = (0.75, 0.05, 0.25, 0.10)
    params_a = (fp = fp, ps = ps, vp = vp, ve = ve)
    params_b = (Tp = vp / fp, Te = ve / ps, vp = vp, ve = ve)
    ct_a = model_exchange(; t, ca, params = params_a)
    ct_b = model_exchange(; t, ca, params = params_b)
    @test ct_a ≈ ct_b
end

@testset "Two compartment filtration model" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    params = (fp = 0.75, ps = 0.05, vp = 0.25, ve = 0.10)
    ct = model_filtration(; t, ca, params)

    est = fit_model(:filtration, :lls; t, ca, ct, mask = [true]).est
    @test round(est.fp[1], digits = 3) == params.fp
    @test round(est.ps[1], digits = 3) == params.ps
    @test round(est.ve[1], digits = 3) == params.ve
    @test round(est.vp[1], digits = 3) == params.vp

    est = fit_model(:filtration, :nls; t, ca, ct, mask = [true]).est
    @test round(est.fp[1], digits = 3) == params.fp
    @test round(est.ps[1], digits = 3) == params.ps
    @test round(est.ve[1], digits = 3) == params.ve
    @test round(est.vp[1], digits = 3) == params.vp

    @test is_all_nan(fit_model(:filtration, :lls; t, ca, ct, mask = false).est)
    @test is_all_nan(fit_model(:filtration, :nls; t, ca, ct, mask = [false]).est)
    @test_throws ErrorException fit_model(:filtration; t, ca, ct, mask = [true, true])

    fp, ps, vp, ve = (0.75, 0.05, 0.25, 0.10)
    params_a = (fp = fp, ps = ps, vp = vp, ve = ve)
    params_b = (Tp = vp / fp, Te = ve / ps, vp = vp, ve = ve)
    ct_a = model_filtration(; t, ca, params = params_a)
    ct_b = model_filtration(; t, ca, params = params_b)
    @test ct_a ≈ ct_b
end

@testset "Reference region models" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    params = (kt = 0.25, kep = 0.5)
    ct = model_tofts(; t, ca, params)
    ref_params = (kt = 0.14, kep = 1.0)
    crr = model_tofts(; t, ca, params = ref_params)

    rel_kt = params.kt / ref_params.kt
    rel_ve = (params.kt / params.kep) / (ref_params.kt / ref_params.kep)
    rel_kt = round(rel_kt, digits = 3)
    rel_ve = round(rel_ve, digits = 3)

    est = fit_model(:rrm, :lls; t, crr, ct, mask = [true]).est
    @test round(est.rel_kt[1], digits = 3) == rel_kt
    @test round(est.rel_ve[1], digits = 3) == rel_ve
    @test round(est.kep[1], digits = 3) == params.kep
    @test round(est.kep_rr[1], digits = 3) == ref_params.kep

    est = fit_model(:rrm, :nls; t, crr, ct).est
    @test round(est.rel_kt[1], digits = 3) == rel_kt
    @test round(est.rel_ve[1], digits = 3) == rel_ve
    @test round(est.kep[1], digits = 3) == params.kep
    @test round(est.kep_rr[1], digits = 3) == ref_params.kep

    est = fit_model(:crrm, :lls; t, crr, ct, mask = [true]).est
    @test round(est.rel_kt[1], digits = 3) == rel_kt
    @test round(est.rel_ve[1], digits = 3) == rel_ve
    @test round(est.kep[1], digits = 3) == params.kep
    @test round(est.kep_rr[1], digits = 3) == ref_params.kep

    est = fit_model(:crrm, :nls; t, crr, ct).est
    @test round(est.rel_kt[1], digits = 3) == rel_kt
    @test round(est.rel_ve[1], digits = 3) == rel_ve
    @test round(est.kep[1], digits = 3) == params.kep
    @test round(est.kep_rr[1], digits = 3) == ref_params.kep

    @test is_all_nan(fit_model(:rrm, :lls; t, crr, ct, mask = false).est)
    @test is_all_nan(fit_model(:rrm, :nls; t, crr, ct, mask = [false]).est)
    @test is_all_nan(fit_model(:crrm, :lls; t, crr, ct, mask = [false]).est)
    @test is_all_nan(fit_model(:crrm, :nls; t, crr, ct, mask = false).est)
    @test_throws ErrorException fit_model(:rrm; t, crr, ct, mask = [true, true])

    # With noise
    num_noisy_replications = 100
    σ = 0.1
    noisy_ct = zeros(num_noisy_replications, length(t))
    for idx = 1:num_noisy_replications
        noisy_ct[idx, :] .= ct .+ σ .* randn(length(ct))
    end

    rel_kt = round(rel_kt, digits = 1)
    rel_ve = round(rel_ve, digits = 1)

    est = fit_model(:rrm, :lls; t, crr, ct = noisy_ct).est
    @test round(mean(est.rel_kt), digits = 1) == rel_kt
    @test round(mean(est.rel_ve), digits = 1) == rel_ve
    @test round(mean(est.kep), digits = 1) == params.kep
    @test round(mean(est.kep_rr), digits = 1) == ref_params.kep

    est_constrained = fit_model(:crrm, :lls; t, crr, ct = noisy_ct).est
    @test round(mean(est_constrained.rel_kt), digits = 1) == rel_kt
    @test round(mean(est_constrained.rel_ve), digits = 1) == rel_ve
    @test round(mean(est_constrained.kep), digits = 1) == params.kep
    @test round(mean(est_constrained.kep_rr), digits = 1) == ref_params.kep
    # Constrained est should have better precision/variability than unconstrained
    @test std(est_constrained.rel_kt) < std(est.rel_kt)
    @test std(est_constrained.rel_ve) < std(est.rel_ve)
    @test std(est_constrained.kep) < std(est.kep)
end

@testset "In-vivo demo" begin
    destination = "./demo_data"
    (vfa_folder, dce_folder) = Perfusion.download_invivo_studies(; destination)
    vfa = load_vfa_dicom(vfa_folder)
    dce = load_dce_dicom(dce_folder; num_slices = 16)

    relaxation_maps = fit_relaxation(:despot; vfa...).est

    bolus_arrival_frame = 3
    r1 = 3.3 / 1000
    concentration = signal_to_concentration(
        dce.signal;
        R10 = 1.0 ./ relaxation_maps.T1,
        angle = dce.angle,
        TR = dce.TR,
        r1 = r1,
        BAF = bolus_arrival_frame,
    )

    Perfusion.download_invivo_studies(; destination) # for more code coverage
end
