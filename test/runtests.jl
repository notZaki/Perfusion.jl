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

@testset "Tofts model" begin
    t = collect(1:600) ./ 60
    Cp = aif_parker(t .- 1)

    test_params = (ktrans = 0.25, kep = 0.5)
    Ct = model_tofts(t=t, Cp=Cp, parameters=test_params)

    estimates = fit_tofts(t=t, Cp=Cp, Ct=Ct, mask=[true]).estimates
    @test round(estimates.ktrans[1], digits=3) == test_params.ktrans
    @test round(estimates.kep[1], digits=3) == test_params.kep

    estimates = fit_tofts(t=t, Cp=Cp, Ct=Ct, method=:nls).estimates
    @test round(estimates.ktrans[1], digits=3) == test_params.ktrans
    @test round(estimates.kep[1], digits=3) == test_params.kep

    @test fit_tofts(t=t, Cp=Cp, Ct=Ct, mask=false).estimates == (ktrans=[0.0], kep=[0.0])
    @test fit_tofts(t=t, Cp=Cp, Ct=Ct, mask=false, method=:nls).estimates == (ktrans=[0.0], kep=[0.0])
    @test_throws ErrorException fit_tofts(t=t, Cp=Cp, Ct=Ct, method=:NLLS)
    @test_throws ErrorException fit_tofts(t=t, Cp=Cp, Ct=Ct, mask=[true, true])
end

@testset "Extended Tofts model" begin
    t = collect(1:600) ./ 60
    Cp = aif_georgiou(t .- 1)

    test_params = (ktrans = 0.5, kep = 1.0, vp = 0.1)
    Ct = model_tofts(t=t, Cp=Cp, parameters=test_params)

    estimates = fit_extendedtofts(t=t, Cp=Cp, Ct=Ct, mask=[true]).estimates
    @test round(estimates.ktrans[1], digits=3) == test_params.ktrans
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.vp[1], digits=3) == test_params.vp

    estimates = fit_extendedtofts(t=t, Cp=Cp, Ct=Ct, method=:nls).estimates
    @test round(estimates.ktrans[1], digits=3) == test_params.ktrans
    @test round(estimates.kep[1], digits=3) == test_params.kep
    @test round(estimates.vp[1], digits=3) == test_params.vp

    @test fit_extendedtofts(t=t, Cp=Cp, Ct=Ct, mask=false).estimates == (ktrans=[0.0], kep=[0.0], vp=[0.0])
    @test fit_extendedtofts(t=t, Cp=Cp, Ct=Ct, mask=false, method=:nls).estimates == (ktrans=[0.0], kep=[0.0], vp=[0.0])
    @test_throws ErrorException fit_extendedtofts(t=t, Cp=Cp, Ct=Ct, method=:NLLS)
    @test_throws ErrorException fit_extendedtofts(t=t, Cp=Cp, Ct=Ct, mask=[true, true])
end
