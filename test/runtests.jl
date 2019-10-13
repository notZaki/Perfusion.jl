using Perfusion
using Test
using UnicodePlots

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
        "Weinmann" => aif_weinmann
    )
    plot_cfg(name) = (title=name, xlabel = "Gd [mM]", ylabel = "time [min]", canvas = DotCanvas)
    printdiv() = printstyled("\n==================================================================================\n", color=:light_black);

    for (name, func) in aifs_to_test
        println( lineplot(scan_timepoints, func(aif_timepoints); plot_cfg(name)...) )
        printdiv()
    end
end

function is_all_nan(x::NamedTuple)
    for val in x
        if any(@. !isnan(val))
            return false
        end
    end
    return true
end

@testset "Tofts model" begin
    t = collect(1:600) ./ 60
    cp = aif_parker(t .- 1)

    test_params = (kt = 0.25, kep = 0.5)
    ct = model_tofts(t=t, cp=cp, parameters=test_params)

    estimates = fit_model(:tofts, :lls, t=t, cp=cp, ct=ct, mask=[true]).estimates
    @test round(estimates.kt[1], digits=3) == test_params.kt
    @test round(estimates.kep[1], digits=3) == test_params.kep

    estimates = fit_model(:tofts, t=t, cp=cp, ct=ct).estimates
    @test round(estimates.kt[1], digits=3) == test_params.kt
    @test round(estimates.kep[1], digits=3) == test_params.kep

    @test is_all_nan(fit_model(:tofts, :lls, t=t, cp=cp, ct=ct, mask=false).estimates)
    @test is_all_nan(fit_model(:tofts, :nls, t=t, cp=cp, ct=ct, mask=[false]).estimates)
    @test_throws ErrorException fit_model(:tofts, t=t, cp=cp, ct=ct, mask=[true, true])
end

@testset "Extended Tofts model" begin
    t = collect(1:600) ./ 60
    cp = aif_georgiou(t .- 1)

    test_params = (kt = 0.5, kep = 1.0, vp = 0.1)
    ct = model_tofts(t=t, cp=cp, parameters=test_params)

    estimates = fit_model(:extendedtofts, :lls, t=t, cp=cp, ct=ct, mask=[true]).estimates
    @test round(estimates.kt[1], digits=3) == test_params.kt
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.vp[1], digits=3) == test_params.vp

    estimates = fit_model(:extendedtofts, t=t, cp=cp, ct=ct).estimates
    @test round(estimates.kt[1], digits=3) == test_params.kt
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.vp[1], digits=3) == test_params.vp

    @test is_all_nan(fit_model(:extendedtofts, :lls, t=t, cp=cp, ct=ct, mask=false).estimates)
    @test is_all_nan(fit_model(:extendedtofts, :nls, t=t, cp=cp, ct=ct, mask=[false]).estimates)
    @test_throws ErrorException fit_model(:extendedtofts, t=t, cp=cp, ct=ct, mask=[true, true])
end

@testset "Compartmental tissue uptake model" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    test_params = (fp = 0.75, ps = 0.05, vp = 0.25)
    ct = model_uptake(t=t, ca=ca, parameters=test_params)

    estimates = fit_model(:uptake, :lls, t=t, ca=ca, ct=ct, mask=[true]).estimates
    @test round(estimates.fp[1], digits=2) == test_params.fp
    @test round(estimates.ps[1], digits=2) == test_params.ps
    @test round(estimates.vp[1], digits=2) == test_params.vp

    estimates = fit_model(:uptake, :nls, t=t, ca=ca, ct=ct, mask=[true]).estimates
    @test round(estimates.fp[1], digits=3) == test_params.fp
    @test round(estimates.ps[1], digits=3) == test_params.ps
    @test round(estimates.vp[1], digits=3) == test_params.vp

    @test is_all_nan(fit_model(:uptake, :lls, t=t, ca=ca, ct=ct, mask=false).estimates)
    @test is_all_nan(fit_model(:uptake, :nls, t=t, ca=ca, ct=ct, mask=[false]).estimates)
    @test_throws ErrorException fit_model(:uptake, t=t, ca=ca, ct=ct, mask=[true, true])
end

@testset "Two compartment exchange model" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    test_params = (fp = 0.75, ps = 0.05, vp = 0.25, ve = 0.10)
    ct = model_exchange(t=t, ca=ca, parameters=test_params)

    estimates = fit_model(:exchange, :lls, t=t, ca=ca, ct=ct, mask=[true]).estimates
    @test round(estimates.fp[1], digits=3) == test_params.fp
    @test round(estimates.ps[1], digits=3) == test_params.ps
    @test round(estimates.ve[1], digits=3) == test_params.ve
    @test round(estimates.vp[1], digits=3) == test_params.vp

    estimates = fit_model(:exchange, :nls, t=t, ca=ca, ct=ct, mask=[true]).estimates
    @test round(estimates.fp[1], digits=3) == test_params.fp
    @test round(estimates.ps[1], digits=3) == test_params.ps
    @test round(estimates.ve[1], digits=3) == test_params.ve
    @test round(estimates.vp[1], digits=3) == test_params.vp

    @test is_all_nan(fit_model(:exchange, :lls, t=t, ca=ca, ct=ct, mask=false).estimates)
    @test is_all_nan(fit_model(:exchange, :nls, t=t, ca=ca, ct=ct, mask=[false]).estimates)
    @test_throws ErrorException fit_model(:exchange, t=t, ca=ca, ct=ct, mask=[true, true])

    fp, ps, vp, ve = (0.75, 0.05, 0.25, 0.10)
    params_a = (fp = fp, ps = ps, vp = vp, ve = ve)
    params_b = (Tp = vp/fp, Te = ve/ps, vp = vp, ve=ve)
    ct_a = model_exchange(t=t, ca=ca, parameters=params_a)
    ct_b = model_exchange(t=t, ca=ca, parameters=params_b)
    @test ct_a ≈ ct_b
end

@testset "Two compartment filtration model" begin
    t = collect(1:600) ./ 60
    ca = aif_georgiou(t .- 1)

    test_params = (fp = 0.75, ps = 0.05, vp = 0.25, ve = 0.10)
    ct = model_filtration(t=t, ca=ca, parameters=test_params)

    estimates = fit_model(:filtration, :lls, t=t, ca=ca, ct=ct, mask=[true]).estimates
    @test round(estimates.fp[1], digits=3) == test_params.fp
    @test round(estimates.ps[1], digits=3) == test_params.ps
    @test round(estimates.ve[1], digits=3) == test_params.ve
    @test round(estimates.vp[1], digits=3) == test_params.vp

    estimates = fit_model(:filtration, :nls, t=t, ca=ca, ct=ct, mask=[true]).estimates
    @test round(estimates.fp[1], digits=3) == test_params.fp
    @test round(estimates.ps[1], digits=3) == test_params.ps
    @test round(estimates.ve[1], digits=3) == test_params.ve
    @test round(estimates.vp[1], digits=3) == test_params.vp

    @test is_all_nan(fit_model(:filtration, :lls, t=t, ca=ca, ct=ct, mask=false).estimates)
    @test is_all_nan(fit_model(:filtration, :nls, t=t, ca=ca, ct=ct, mask=[false]).estimates)
    @test_throws ErrorException fit_model(:filtration, t=t, ca=ca, ct=ct, mask=[true, true])

    fp, ps, vp, ve = (0.75, 0.05, 0.25, 0.10)
    params_a = (fp = fp, ps = ps, vp = vp, ve = ve)
    params_b = (Tp = vp/fp, Te = ve/ps, vp = vp, ve=ve)
    ct_a = model_filtration(t=t, ca=ca, parameters=params_a)
    ct_b = model_filtration(t=t, ca=ca, parameters=params_b)
    @test ct_a ≈ ct_b
end

@testset "Reference region models" begin
    t = collect(1:600) ./ 60
    cp = aif_georgiou(t .- 1)

    test_params = (kt = 0.25, kep = 0.5)
    ct = model_tofts(t=t, cp=cp, parameters=test_params)
    ref_params = (kt = 0.14, kep = 1.0)
    crr = model_tofts(t=t, cp=cp, parameters=ref_params)

    rel_kt = test_params.kt / ref_params.kt
    rel_ve = (test_params.kt / test_params.kep) / (ref_params.kt / ref_params.kep)
    rel_kt = round(rel_kt, digits=3)
    rel_ve = round(rel_ve, digits=3)

    estimates = fit_model(:referenceregion, :lls, t=t, crr=crr, ct=ct, mask=[true]).estimates
    @test round(estimates.rel_kt[1], digits=3) == rel_kt
    @test round(estimates.rel_ve[1], digits=3) == rel_ve
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.kep_rr[1], digits=3) == ref_params.kep

    estimates = fit_model(:referenceregion, :nls, t=t, crr=crr, ct=ct).estimates
    @test round(estimates.rel_kt[1], digits=3) == rel_kt
    @test round(estimates.rel_ve[1], digits=3) == rel_ve
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.kep_rr[1], digits=3) == ref_params.kep

    estimates = fit_model(:constrained_referenceregion, :lls, t=t, crr=crr, ct=ct, mask=[true]).estimates
    @test round(estimates.rel_kt[1], digits=3) == rel_kt
    @test round(estimates.rel_ve[1], digits=3) == rel_ve
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.kep_rr[1], digits=3) == ref_params.kep

    estimates = fit_model(:constrained_referenceregion, :nls, t=t, crr=crr, ct=ct).estimates
    @test round(estimates.rel_kt[1], digits=3) == rel_kt
    @test round(estimates.rel_ve[1], digits=3) == rel_ve
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.kep_rr[1], digits=3) == ref_params.kep

    @test is_all_nan(fit_model(:referenceregion, :lls, t=t, crr=crr, ct=ct, mask=false).estimates)
    @test is_all_nan(fit_model(:referenceregion, :nls, t=t, crr=crr, ct=ct, mask=[false]).estimates)
    @test is_all_nan(fit_model(:constrained_referenceregion, :lls, t=t, crr=crr, ct=ct, mask=[false]).estimates)
    @test is_all_nan(fit_model(:constrained_referenceregion, :nls, t=t, crr=crr, ct=ct, mask=false).estimates)
    @test_throws ErrorException fit_model(:referenceregion, t=t, crr=crr, ct=ct, mask=[true, true])
end
