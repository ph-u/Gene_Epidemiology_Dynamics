# func.jl — shared utility functions
# Loaded by: nfds.jl  model.jl  tests/runtests.jl
# Standard library only; no external packages needed here.

# ─── Configuration reader ──────────────────────────────────────────────────────

"""
Read one value from input.csv by its "Type" key.
Returns `nothing` when the key is absent.

    read_cfg("../data/input.csv", "population")  # → "800000"
"""
function read_cfg(csv_path::String, key::String)
    for (i, line) in enumerate(eachline(csv_path))
        i == 1 && continue                           # skip header
        cols = split(line, ','; limit = 3)
        length(cols) >= 2 && strip(cols[1]) == key && return strip(cols[2])
    end
    return nothing
end

# ─── Jensen-Shannon Divergence ─────────────────────────────────────────────────

"""
Jensen-Shannon divergence between two count vectors p and q.

Analogy: p and q are histograms (e.g. which genotypes were seen at a time-point).
JSD measures how different their shapes are.

  JSD = 0.0       → identical shapes
  JSD = log(2)    → no overlap at all   (≈ 0.693)
  JSD = NaN       → one vector is all zeros — no data to compare

Counts are normalised internally; total sizes do not matter.
"""
function jsd(p::AbstractVector{<:Real}, q::AbstractVector{<:Real})
    sp = sum(p);  sq = sum(q)
    (sp == 0 || sq == 0) && return NaN
    p̂ = p ./ sp;  q̂ = q ./ sq          # convert counts → proportions (sum to 1)
    m = (p̂ .+ q̂) ./ 2                   # midpoint distribution
    # KL divergence: extra bits to describe distribution a using code for m
    kl(a, b) = sum(aᵢ > 0 ? aᵢ * log(aᵢ / bᵢ) : 0.0 for (aᵢ, bᵢ) in zip(a, b))
    return 0.5kl(p̂, m) + 0.5kl(q̂, m)
end

"""
Sum JSD over all sampling time-points (columns of two equal-shaped matrices).
Empty columns are penalised by log(2), the maximum possible distance.
"""
function total_jsd(sim::AbstractMatrix, obs::AbstractMatrix)
    @assert size(sim) == size(obs)  "Shape mismatch: $(size(sim)) ≠ $(size(obs))"
    return sum(let d = jsd(view(sim, :, t), view(obs, :, t))
                   isnan(d) ? log(2) : d
               end for t in 1:size(obs, 2))
end

# ─── Weighted sampling ─────────────────────────────────────────────────────────

"""
Draw `n` indices from 1:length(weights) with probability proportional to weights.
Sampling is with replacement.  Weights need not sum to 1.
Equivalent to R's  sample(1:k, n, replace=TRUE, prob=weights).
"""
function wt_sample(weights::AbstractVector{<:Real}, n::Int)
    cdf = cumsum(weights)
    cdf ./= cdf[end]           # normalise (guards against floating-point rounding)
    return [searchsortedfirst(cdf, rand()) for _ in 1:n]
end
