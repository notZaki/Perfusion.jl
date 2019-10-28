
function _fitRRM(;t::AbstractVector, Ct::AbstractArray, Crr::AbstractVector)
    initial_size_Ct = size(Ct)[1:end-1]
    Ct = _melt_array_to_2D(Ct)
    (numVox, sT) = size(Ct)

    pkParams = zeros(numVox,3)
    resid = zeros(numVox)

    M = zeros(sT,3)
    M[:,1] = Crr
    M[:,2] = cumtrapz(t, Crr)

    for i=1:numVox
            M[:,3] .= -cumtrapz(t, Ct[i,:])
            pkParams[i,:] = M\Ct[i,:]
            resid[i] = Statistics.mean((M*pkParams[i,:]-Ct[i,:]).^2)
    end

    pkParams = _restore_array(pkParams, initial_size_Ct)
    resid = _restore_array(resid, initial_size_Ct)
    return(params=pkParams, resid=resid)
end

function fitCRRM(;t::AbstractVector, Ct::AbstractArray, Crr::AbstractVector)
    initial_size_Ct = size(Ct)[1:end-1]
    Ct = _melt_array_to_2D(Ct)
    (numVox, sT) = size(Ct)


    estRRM = _fitRRM(t=t, Ct=Ct, Crr=Crr)
    kepRR = _estimate_kepRR(estRRM.params)


    pkParams = zeros(numVox,3)
    resid = zeros(numVox)

    M = zeros(sT,2)
    M[:,1] = Crr + kepRR * cumtrapz(t, Crr)
    myCatch = []
    for i=1:numVox
        M[:,2] .= -cumtrapz(t, Ct[i,:])
        pkParams[i,1:2] = M\Ct[i,:]
        resid[i] = Statistics.mean((M*pkParams[i,1:2]-Ct[i,:]).^2)
    end
    pkParams[:,3] .= pkParams[:,2]
    pkParams[:,2] .= kepRR * pkParams[:,1] ./ pkParams[:,3]


    pkParams = _restore_array(pkParams, initial_size_Ct)
    resid = _restore_array(resid, initial_size_Ct)
    return(params=pkParams, resid=resid, kepRR=kepRR, myCatch=myCatch, estRRM=estRRM)
end

function RRIFT(; Cp::AbstractVector, Crr::AbstractVector, t::AbstractVector, kepRR::Number, tail_start::Int=0)
    if tail_start == 0
        tail_start = 1
        tail_frames = 1:length(t)
    else
        tail_frames = tail_start:length(t)
    end
    num = Crr[tail_frames] .- Crr[tail_start] .+ kepRR .* cumtrapz(t[tail_frames], Crr[tail_frames])
    denum = cumtrapz(t[tail_frames], Cp[tail_frames])
    return denum\num
end



function _estimate_kepRR(rrm_params)
    goodVals = vec( all(rrm_params .> 0, dims = 2) )
    interquartile_mean( rrm_params[goodVals,2] ./ rrm_params[goodVals,1] )
end



function get_self_reference_region(Ct::AbstractArray, num_clusters::Int=3)
    candidates = nmf_hierarchical(_melt_array_to_2D(Ct)', num_clusters)
    num_timepoints = size(candidates)[1]
    num_initial_timepoints = Int(round(num_timepoints/3))
    initial_timepoints = collect(1:num_initial_timepoints)
    initial_area_under_curve = zeros(num_clusters)
    for i=1:num_clusters
        initial_area_under_curve[i] = trapz(initial_timepoints, candidates[initial_timepoints,i])
    end
    return candidates[:, sortperm(initial_area_under_curve)]
end
