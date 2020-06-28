function expconv(A::AbstractVector, B::Number, t::AbstractVector)
    # Returns f = A ConvolvedWith exp(-B t)
    # Based on Flouri et al. (2016) MRM 76(3), doi: 10.1002/mrm.25991
    @assert length(A) == length(t)
    f = zeros(length(t))
    for i = 2:length(t)
        x = B * (t[i] - t[i-1])
        dA = (A[i] - A[i-1]) / x
        E = exp(-x)
        E0 = 1 - E
        E1 = x - E0
        f[i] = E * f[i-1] + A[i-1] * E0 + dA * E1
    end
    return f ./ B
end

macro extract(varnames, namedtuple)
    ex = Expr(:block)
    ex.args = [:($(esc(var)) = getindex($(esc(namedtuple)), $(esc(QuoteNode(var))))) for var in varnames.args]
    ex
end

function interquartile_mean(x::AbstractVector)
    if std(x) > 1e-3
        quartiles = quantile(x, [0.25, 0.75])
        # `<=` is safer than `<` because latter can creat empty array if x is short
        interquartile_x = @. x[quartiles[1] <= x <= quartiles[2]]
    else
        interquartile_x = x
    end
    return mean(interquartile_x)
end

function percent_error(estimate::Number, truth::Number)
    return 100 * (estimate - truth) / truth
end

function positive_only_mask(x::NamedTuple)
    valid_keys = keys(x)
    mask = trues(size(x[valid_keys[1]]))
    for array in x
        @. mask = mask & (array > 0)
    end
    return mask
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

function unify_size(element, desired_size)
    if size(element) == desired_size
        return element
    elseif size(element) == ()
        return repeat([element], desired_size...)
    else
        error("Element size: $(size(element)) does not match input size $(desired_size)")
    end
end

