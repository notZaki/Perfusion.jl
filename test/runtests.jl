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
