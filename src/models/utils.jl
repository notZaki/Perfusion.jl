function interquartile_mean(x)
    x = vec(x)
    if std(x)>1e-3
        quartiles = quantile(x,[0.25, 0.75])
        x = x[x.>quartiles[1]]
        x = x[x.<quartiles[2]]
    end
    x = x[.!isnan.(x)]
    return mean(x)
end

function percenterror(estV::Number, trueV::Number)::Float64
    return 100 * (estV - trueV) / trueV
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
