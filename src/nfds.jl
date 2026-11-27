# nfds.jl — NFDS population genetics model
# Equivalent to R's setup.r + nfds.r combined.
#
# Two entry points used by model.jl:
#   build_state(input_csv)  → ModelState   (call ONCE — heavy I/O + computation)
#   abc_distance(θ, state)  → Float64      (call MANY TIMES in the ABC loop)
#
# ── PLAIN-ENGLISH MODEL DESCRIPTION ──────────────────────────────────────────
# We simulate ~800 000 bacteria living inside a human host population.
# Every bacterium belongs to a "genotype" — a unique 0/1 pattern across ~400 genes.
# Each simulated month:
#   A. Count how many of each genotype currently exist (in each sub-population).
#   B. Compute NFDS fitness: genotypes that are RARER than their historical average
#      are favoured → diversity is maintained over time.  (KEY matrix operation ★)
#   C. Apply vaccine pressure: bacteria targeted by the vaccine reproduce less.
#   D. Sample offspring (Poisson or Negative-Binomial stochastic reproduction).
#   E. Allow bacteria to migrate between sub-populations.
# At sampling months we record genotype counts and compare to real survey data
# using Jensen-Shannon Divergence (JSD).  The ABC algorithm searches for the
# parameter values that minimise this distance.

# @isdefined guard: if func.jl was already included (e.g. from runtests.jl),
# skip re-loading it — redefining structs is a hard error in Julia.
# joinpath(@__DIR__, ...) resolves relative to this file's directory (src/),
# not the working directory, so the path is correct from any call site.
@isdefined(jsd) || include(joinpath(@__DIR__, "func.jl"))

using CSV, DataFrames, Distributions, StatsBase, JLD2
using LinearAlgebra: mul!
using Statistics: mean
using Random: rand, randn, randperm, seed!
using Printf: @printf

# ─── Data structure ────────────────────────────────────────────────────────────

"""
ModelState: all pre-computed quantities.
Build once with build_state(); pass to abc_distance() as many times as needed.
Fields are documented with their shapes; nG=genotypes, nL=loci,
nD=demes, nV=vaccines, nT=sampling time-points, nA=age levels.
"""
struct ModelState
    # ── Genotype lookup tables ─────────────────────────────────────────────
    # G0[g, l] = 1 if genotype g carries gene l, else 0.
    # Think of G0 as a dictionary: "which genes does this genotype have?"
    G0      ::Matrix{Float64}   # (nG × nL) binary
    VTmat   ::Matrix{Float64}   # (nG × nV) binary — 1 = vaccine-targeted type

    # ── Observed data (what we try to match) ───────────────────────────────
    mIg0    ::Matrix{Int}       # (nG*nD × nT) observed genotype counts
    nObs    ::Matrix{Int}       # (nT × nD)   sample sizes per timepoint per deme
    eqm_dates::Vector{Int}     # which months were sampled, e.g. [36, 72]
    n_gen   ::Int               # last sampling month = simulation length

    # ── Pre-vaccine equilibrium frequencies ────────────────────────────────
    eqm_pre ::Matrix{Float64}  # (nL × nD) target gene freq per deme before vaccine
    eqm_glb ::Vector{Float64}  # (nL,)     global target gene freq
    sel_rank ::Vector{Int}     # gene indices sorted by selection strength (rank 1 = strongest)

    # ── Population structure ───────────────────────────────────────────────
    kD      ::Vector{Int}       # (nD,) host population per deme
    dProp   ::Vector{Float64}  # (nD,) proportion of isolates per deme (sums to 1)
    mP0     ::Matrix{Float64}  # (nG × nD) migrant pool, SC-balanced per deme
    tag_pre ::Vector{Vector{Int}}  # nD × (varies) pre-vaccine genotype pool

    # ── Vaccine rollout schedule ───────────────────────────────────────────
    vstart  ::Matrix{Float64}  # (nV × nD) start month; Inf = never introduced
    a_coh   ::Matrix{Float64}  # (nV × nD) max host age (months) for vaccination
    u_pt    ::Matrix{Float64}  # (nV × nD) uptake fraction 0–1

    # ── Host age structure ─────────────────────────────────────────────────
    a_lev   ::Vector{Int}       # possible host ages in months, e.g. [6,18,36,60]
    a_prob  ::Matrix{Float64}  # (nA × nD) probability of each age in each deme

    # ── Dimensions ────────────────────────────────────────────────────────
    nD ::Int; nV ::Int; nL ::Int; nG ::Int
end

# ─── State build (equiv. setup.r) ─────────────────────────────────────────────

"""
Read and pre-process all data files.
Call ONCE and pass the result to abc_distance() many times.

`input_csv` = path to data/input.csv (the master config).
"""
function build_state(input_csv::String)
    g(key) = read_cfg(input_csv, key)   # shorthand: g("population") → "800000"

    println("Reading data files…")
    d0   = CSV.read(g("data"),    DataFrame; delim = '\t')
    roll = CSV.read(g("rollout"), DataFrame; delim = '\t')
    k    = parse(Int,     g("population"))
    _p0  = parse(Float64, g("percentage initial infected")) / 100

    # ── Identify columns ──────────────────────────────────────────────────
    sc_pos    = findfirst(==("SC"), names(d0))
    vt_cols   = filter(n -> occursin(r"^VT\d+$", n), names(d0))
    gene_cols = names(d0)[sc_pos+1:end]
    deme_lvls = sort(unique(d0.Deme))
    nD = length(deme_lvls);  nV = length(vt_cols);  nL = length(gene_cols)
    d0.deme_id = [findfirst(==(x), deme_lvls) for x in d0.Deme]

    # ── Unique genotypes ──────────────────────────────────────────────────
    # Two isolates with the same 0/1 gene pattern are the same genotype.
    d0.tag     = [join(Vector(row[gene_cols]), "") for row in eachrow(d0)]
    tags_u     = unique(d0.tag)
    nG         = length(tags_u)
    d0.geno_id = [findfirst(==(t), tags_u) for t in d0.tag]
    println("  $nG unique genotypes | $nL loci | $nD demes | $nV vaccines")

    # Build G0: for each unique genotype, pick its first occurrence to get the gene row.
    # G0[g, l] = 1 if genotype g has gene l.
    first_occ = [findfirst(==(i), d0.geno_id) for i in 1:nG]
    G0    = Float64[d0[first_occ[g], c] for g in 1:nG, c in gene_cols]
    VTmat = Float64[d0[first_occ[g], v] for g in 1:nG, v in vt_cols]

    # ── Pre-vaccine equilibrium ───────────────────────────────────────────
    pre = d0[d0.Time .<= 0, :]
    # eqm_pre[l, d] = average frequency of gene l in pre-vaccine isolates from deme d.
    # Fallback to global mean when a deme has no pre-vaccine isolates.
    eqm_glb = Float64[mean(pre[:, l]) for l in gene_cols]
    eqm_pre = Float64[
        let rows = filter(r -> r.deme_id == d, pre)
            isempty(rows) ? eqm_glb[li] : mean(rows[:, l])
        end
        for (li, l) in enumerate(gene_cols), d in 1:nD]

    # Rank genes by how much their global frequency changed after vaccine introduction.
    # High change → strong NFDS signal → should be under stronger selection.
    post      = d0[d0.Time .> 0, :]
    post_freq = Float64[isempty(post) ? eqm_glb[li] : mean(post[:, l])
                        for (li, l) in enumerate(gene_cols)]
    sel_rank  = sortperm(abs.(eqm_glb .- post_freq); rev = true)

    # ── Sampling timepoints and observed count matrix ─────────────────────
    eqm_dates = sort(filter(>(0), unique(d0.Time)))
    nT        = length(eqm_dates)
    n_gen     = maximum(eqm_dates)

    # mIg0 rows: (deme-1)*nG + genotype_index.  Columns: sampling timepoints.
    mIg0 = zeros(Int, nG * nD, nT)
    nObs = zeros(Int, nT, nD)
    for (ti, t) in enumerate(eqm_dates)
        for row in eachrow(filter(r -> r.Time == t, d0))
            g_id, d_id = row.geno_id, row.deme_id
            mIg0[(d_id-1)*nG + g_id, ti] += 1
            nObs[ti, d_id] += 1
        end
    end

    # ── Population sizes ──────────────────────────────────────────────────
    dProp = Float64[count(==(d), d0.deme_id) / nrow(d0) for d in 1:nD]
    kD    = max.(1, round.(Int, k .* dProp))

    # ── Migrant pool (SC-balanced per deme) ───────────────────────────────
    # External migrants should represent all SC labels equally, not just common ones.
    mP0 = zeros(Float64, nG, nD)
    for d in 1:nD
        d_rows  = filter(r -> r.deme_id == d, d0)
        sc_uniq = unique(d_rows.SC)
        for sc in sc_uniq
            sc_genos = unique(filter(r -> r.SC == sc, d_rows).geno_id)
            w = 1.0 / (length(sc_uniq) * length(sc_genos))
            for gi in sc_genos; mP0[gi, d] += w; end
        end
        # Normalise to guard against floating-point drift
        s = sum(mP0[:, d]);  s > 0 && (mP0[:, d] ./= s)
    end

    # ── Vaccine rollout schedule ──────────────────────────────────────────
    vstart = fill(Inf, nV, nD)
    a_coh  = zeros(Float64, nV, nD)
    u_pt   = zeros(Float64, nV, nD)
    for row in eachrow(roll)
        vi = findfirst(==(string(row.Vaccine)), vt_cols)
        di = findfirst(==(string(row.Deme)),    deme_lvls)
        (isnothing(vi) || isnothing(di)) && continue
        vstart[vi, di] = Float64(row.StartMonth)
        a_coh[vi, di]  = Float64(row.AgeWindow)
        u_pt[vi, di]   = Float64(row.Uptake)
    end

    # ── Host age distribution ─────────────────────────────────────────────
    a_lev_raw = sort(unique(d0.Age))
    a_lev     = convert(Vector{Int}, a_lev_raw)
    nA        = length(a_lev)
    a_prob = Float64[
        let n_tot = count(==(d), d0.deme_id)
            n_tot == 0 ? 1.0/nA : count((d0.deme_id .== d) .& (d0.Age .== a)) / n_tot
        end
        for a in a_lev, d in 1:nD]

    # ── Pre-vaccine isolation pool ────────────────────────────────────────
    tag_pre = Vector{Vector{Int}}(undef, nD)
    for d in 1:nD
        pool = d0[(d0.Time .<= 0) .& (d0.deme_id .== d), :geno_id]
        tag_pre[d] = isempty(pool) ? collect(1:nG) : pool  # fallback if no pre-vac data
    end

    println("  Sparsity: $(round(100*mean(==(0), mIg0); digits=1))% zero cells")
    return ModelState(G0, VTmat, mIg0, nObs, eqm_dates, n_gen,
                      eqm_pre, eqm_glb, sel_rank,
                      kD, dProp, mP0, tag_pre,
                      vstart, a_coh, u_pt,
                      a_lev, a_prob,
                      nD, nV, nL, nG)
end

"""Save a ModelState to disk (JLD2 / HDF5 format, loads back in <1 s)."""
save_state(path::String, s::ModelState) = jldsave(path; state = s)
load_state(path::String)                = jldopen(f -> f["state"], path)

# ─── Forward simulation (equiv. m.nfds in nfds.r) ────────────────────────────

"""
Simulate one forward run of the NFDS model for parameter vector θ.

θ is a vector in [0,1]^d (all parameters live in the unit interval and are
rescaled to biological units internally — see model.jl for the layout).

Returns the simulated genotype-count matrix (same shape as state.mIg0),
or `nothing` if the population went extinct or exploded (invalid run).
"""
function simulate(θ::AbstractVector{<:Real}, s::ModelState)
    nD = s.nD;  nV = s.nV;  nL = s.nL;  nG = s.nG

    # ── Rescale parameters to biological units ────────────────────────────
    prop_strong = θ[1]            # fraction of loci under strong NFDS  (0 → 1)
    f_sel       = θ[2] * 0.22    # strong-locus NFDS strength           (0 → 0.22/month)
    w_sel       = θ[3] * 0.01    # weak-locus  NFDS strength            (0 → 0.01/month)
    migration   = θ[4] * 0.2     # external immigration rate            (0 → 0.2 /month)
    co_inf      = θ[5]            # 0 = Poisson offspring; >0 = NegBin overdispersion

    nMW = nD * (nD - 1)
    mW  = zeros(Float64, nD, nD)   # within-deme migration matrix; diagonal is always 0
    pair = 1
    for d1 in 1:nD, d2 in 1:nD
        d1 == d2 && continue
        mW[d1, d2] = θ[5 + pair] * 0.1   # max 0.1 / month per directed pair
        pair += 1
    end
    v_sel = [θ[5 + nMW + v] * 0.5 for v in 1:nV]   # max 0.5 / month
    lVs   = log1p.(-v_sel)                            # log-fitness reduction per vaccine

    # ── Selection strength per gene ────────────────────────────────────────
    # The top `prop_strong` fraction of genes (by historical frequency change,
    # pre-ranked in sel_rank) get the stronger coefficient f_sel; the rest get w_sel.
    n_strong = max(1, round(Int, prop_strong * nL))
    sEl      = fill(w_sel, nL)
    sEl[s.sel_rank[1:n_strong]] .= f_sel

    # ── Initialise population ─────────────────────────────────────────────
    # gen[i] = genotype index of bacterium i (1..nG)
    # dem[i] = deme index     of bacterium i (1..nD)
    # We initialise by drawing kD[d] bacteria from the pre-vaccine pool of deme d.
    gen = Int[];  dem = Int[]
    for d in 1:nD
        isempty(s.tag_pre[d]) && continue
        idx = sample(s.tag_pre[d], s.kD[d]; replace = true)
        append!(gen, idx);  append!(dem, fill(d, s.kD[d]))
    end

    result = zeros(Int, nG * nD, length(s.eqm_dates))
    t_ptr  = 1   # pointer into eqm_dates

    # ═══════════════════════════════════════════════════════════════════════
    for t in 1:s.n_gen   # one loop iteration = one month of bacterial evolution

        n = length(gen)
        n == 0 && return nothing    # population extinct — invalid run

        # ── A: Count how many of each genotype exist in each deme ─────────
        # cnt[g, d] = current number of bacteria of genotype g in deme d.
        cnt  = zeros(Int, nG, nD)
        for i in 1:n;  cnt[gen[i], dem[i]] += 1;  end
        nPer = vec(sum(cnt; dims = 1))   # total bacteria per deme
        any(==(0), nPer) && return nothing

        # ★ B: Compute NFDS fitness (the core matrix operation) ────────────
        #
        # STEP 1 — current gene frequencies:
        #   We need freq[d, l] = fraction of bacteria in deme d that carry gene l.
        #
        #   cnt'    has shape (nD × nG) — transposing to get demes as rows
        #   G0      has shape (nG × nL)
        #   cnt'× G0 is  (nD × nL) — gene COUNTS  per deme (each gene across demes)
        #   Dividing row d by nPer[d] gives gene FREQUENCIES per deme.
        #
        freq = (cnt' * s.G0) ./ reshape(nPer, :, 1)   # (nD × nL)
        #                        ↑ reshape(nPer, nD, 1) broadcasts each row
        #                          separately: row d ÷ nPer[d]

        # STEP 2 — deviation from equilibrium:
        #   dev[d, l] = eqm_pre[l, d] - freq[d, l]
        #   Positive deviation → gene is RARER than target → fitness BOOST
        #   Negative deviation → gene is MORE COMMON than target → fitness PENALTY
        dev = s.eqm_pre' .- freq     # (nD × nL)  s.eqm_pre' is nD×nL

        # STEP 3 — weight each gene's contribution by its selection coefficient:
        #   sEl is (nL,), sEl' is (1 × nL), broadcasts with (nD × nL)
        adj = dev .* sEl'            # (nD × nL)  element-wise, row-broadcast

        # STEP 4 — sum over all genes each genotype carries and exponentiate:
        #   G0 × adj' maps (nG×nL) × (nL×nD) → (nG × nD) fitness per genotype-deme pair.
        #   Taking exp converts additive log-fitness to multiplicative fitness.
        nfds_fit = exp.(s.G0 * adj')   # (nG × nD)

        # For each bacterium, look up its genotype's fitness in its current deme.
        fit = Float64[nfds_fit[gen[i], dem[i]] for i in 1:n]

        # ── C: Assign host ages ────────────────────────────────────────────
        # Every new colonisation event gives the bacterium a fresh host.
        # We draw that host's age from the observed distribution (a_prob).
        h_age = Vector{Int}(undef, n)
        for d in 1:nD
            idx_d  = findall(==(d), dem)
            isempty(idx_d) && continue
            # wt_sample returns indices into a_lev; then we look up the actual age.
            h_age[idx_d] = s.a_lev[wt_sample(s.a_prob[:, d], length(idx_d))]
        end

        # ── D: Vaccine coverage and fitness reduction ──────────────────────
        # Host is covered by vaccine v in deme d if:
        #   (i)  t ≥ vstart[v, d]  — vaccine already introduced, AND
        #   (ii) age at introduction ≤ a_coh[v, d] — host was young enough
        # VT bacteria in covered hosts get a fitness penalty of exp(lVs[v]).
        lvp = zeros(Float64, n)   # log-fitness adjustment per bacterium
        for v in 1:nV, i in 1:n
            d  = dem[i]
            sv = s.vstart[v, d]
            sv == Inf  && continue            # vaccine never came to this deme
            t < sv     && continue            # vaccine not yet started
            h_age[i] - (t - sv) > s.a_coh[v, d] && continue   # host too old
            lvp[i] += s.VTmat[gen[i], v] * lVs[v] * s.u_pt[v, d]
        end
        fit .*= exp.(lvp)

        # ── E: Offspring sampling ──────────────────────────────────────────
        # Each bacterium produces offspring proportional to its fitness.
        # Density scaling: expected total per deme stays near kD[d].
        # coInf = 0 → Poisson (variance = mean, standard).
        # coInf > 0 → NegBin  (variance > mean, more bursty colonisation).
        λ = fit .* Float64[s.kD[dem[i]] / nPer[dem[i]] for i in 1:n]
        offspring = co_inf < 1e-5 ?
            rand.(Poisson.(λ)) :
            [rand(NegativeBinomial(1/co_inf, (1/co_inf)/(1/co_inf + λ[i]))) for i in 1:n]

        total_off = sum(offspring)
        total_off == 0              && return nothing   # extinction
        total_off > 20 * sum(s.kD) && return nothing   # explosion

        # Expand: bacterium i → offspring[i] identical copies.
        # (Pre-allocate the target vectors for speed — avoids repeated re-allocation.)
        new_gen = Vector{Int}(undef, total_off)
        new_dem = Vector{Int}(undef, total_off)
        ptr = 1
        for i in 1:n
            offspring[i] == 0 && continue
            new_gen[ptr:ptr+offspring[i]-1] .= gen[i]
            new_dem[ptr:ptr+offspring[i]-1] .= dem[i]
            ptr += offspring[i]
        end
        gen = new_gen;  dem = new_dem

        # ── F: Within-deme (metapopulation) migration ─────────────────────
        # Host-level movement: when a person moves, all their bacteria move with them.
        # We iterate in random order to avoid directional bias.
        if any(>(0), mW) && nD > 1
            for i in randperm(length(gen))
                d_from = dem[i]
                r_tot  = sum(mW[d_from, :])
                r_tot == 0 && continue
                rand() >= r_tot && continue       # most bacteria do NOT move
                dem[i] = wt_sample(mW[d_from, :] ./ r_tot, 1)[1]
            end
        end

        # ── G: External immigration ────────────────────────────────────────
        # A small fraction of colonisations come from outside the study area.
        # Immigrants are drawn from mP0, which is balanced across SC labels.
        if migration > 0
            n_mig = rand(Poisson(migration * length(gen)))
            for _ in 1:n_mig
                d_to  = wt_sample(s.dProp, 1)[1]
                g_mig = wt_sample(s.mP0[:, d_to], 1)[1]
                push!(gen, g_mig);  push!(dem, d_to)
            end
        end

        # ── H: Record at sampling time-points ─────────────────────────────
        if t_ptr <= length(s.eqm_dates) && t == s.eqm_dates[t_ptr]
            for i in 1:length(gen)
                result[(dem[i]-1)*nG + gen[i], t_ptr] += 1
            end
            t_ptr += 1
        end

    end   # ═══════════════ end generation loop ═══════════════════════════

    return result
end

# ─── ABC distance function ────────────────────────────────────────────────────

"""
Wrapper: run simulate(θ, state) → compute JSD vs observed data.
Called N×20 times per ABC-SMC generation, entirely in memory via Threads.@threads.
Returns the worst possible distance when the simulation fails.
"""
function abc_distance(θ::AbstractVector{<:Real}, s::ModelState)
    sim = simulate(θ, s)
    isnothing(sim) && return Float64(length(s.eqm_dates)) * log(2)
    return total_jsd(Float64.(sim), Float64.(s.mIg0))
end

# ─── Thread-local simulation buffer ──────────────────────────────────────────
#
# Pre-allocating these arrays eliminates the dominant GC pressure source.
# Without buffers, every simulate() call allocates ~100–200 MB of temporaries
# (per-month cnt/freq/dev/adj/nfds_fit matrices, per-bacterium vectors, and
# the result matrix).  With 7 threads × 10,000 proposals × 120 months that is
# ~100–200 GB of cumulative allocation throughput per SMC generation, forcing
# constant GC collections that consume 20–30% of CPU and cause >4 GB heap
# spikes.  With SimBuffer, each of these allocations becomes an in-place write
# into permanently resident memory; GC pressure drops by ~100×.
#
# One SimBuffer per thread.  Call alloc_sim_buffer(state) for each thread and
# store in a Vector{SimBuffer} indexed by Threads.threadid().

mutable struct SimBuffer
    # Fixed-shape matrices — same shape every call, cleared with fill! each use.
    cnt      ::Matrix{Int}       # nG × nD  — genotype counts per deme
    freq     ::Matrix{Float64}   # nD × nL  — gene frequencies
    dev      ::Matrix{Float64}   # nD × nL  — deviation from equilibrium
    adj      ::Matrix{Float64}   # nD × nL  — selection-weighted deviation
    nfds_fit ::Matrix{Float64}   # nG × nD  — NFDS fitness per genotype-deme
    result   ::Matrix{Int}       # nG*nD × n_eqm_dates — simulation output
    result_f ::Matrix{Float64}   # nG*nD × n_eqm_dates — Float64 copy for JSD
    mIg0_f   ::Matrix{Float64}   # nG*nD × n_eqm_dates — Float64(s.mIg0), constant
    sEl      ::Vector{Float64}   # nL — per-gene selection coefficients

    # Pre-computed CDFs (constant across simulate! calls; built once in alloc).
    age_cdf  ::Matrix{Float64}   # nA × nD — CDF for host-age sampling per deme
    dProp_cdf::Vector{Float64}   # nD  — CDF for immigration deme selection
    mP0_cdf  ::Matrix{Float64}   # nG × nD — CDF for migrant genotype per deme

    # Population vectors — pre-allocated at max capacity (20 × sum(kD)).
    # Actual live population is gen[1:n] and dem[1:n]; n tracked as a local.
    # Pointer-swap pattern: after offspring expansion, swap gen↔new_gen and
    # dem↔new_dem (zero-copy rename) rather than allocating fresh arrays.
    gen      ::Vector{Int}
    dem      ::Vector{Int}
    new_gen  ::Vector{Int}
    new_dem  ::Vector{Int}

    # Per-bacterium working vectors (capacity = max_n; sliced to 1:n each month).
    fit      ::Vector{Float64}
    h_age    ::Vector{Int}
    lvp      ::Vector{Float64}
    λ        ::Vector{Float64}
    offspring::Vector{Int}

    # Per-deme index arrays: built in a single O(n) pass each month to avoid
    # repeated findall(==(d), dem_v) which allocates O(n) per deme per month.
    idx_d    ::Vector{Vector{Int}}   # nD arrays, each capacity max_n
    n_per_d  ::Vector{Int}           # nD — live count of bacteria in each deme

    # Permutation buffer for within-deme migration shuffle.
    # randperm(n) allocates n×8 bytes (≈6.4 MB when n≈800,000) every month.
    # With 120 months × 7 threads, that is 5.4 GB of allocation throughput per
    # SMC generation.  The enlarged GC trigger threshold caused by the SimBuffer
    # live set (7×272 MB = 1.9 GB → GC target ≈ 7.8 GB) delays collection,
    # allowing several GB of dead randperm arrays to accumulate — the 7 GB spike.
    # Pre-allocating here and using Random.shuffle!(view) eliminates it entirely.
    perm     ::Vector{Int}           # capacity max_n — reused each month in-place
end

"""
Allocate one SimBuffer sized for the given ModelState.
Call once per thread before entering the ABC-SMC loop.
"""
function alloc_sim_buffer(s::ModelState)
    # Buffer capacity: 5× carrying capacity covers realistic population variation.
    # Simulations that grow beyond this are treated as invalid explosions and
    # return false before writing to the buffer (see bounds checks in simulate!).
    # Using the biology explosion threshold (20×) instead would pre-allocate
    # ~8.5 GB of permanent buffer memory — 7× more than needed — leaving
    # insufficient headroom within the 12 GB LSF memory limit.
    max_n = 5 * sum(s.kD)
    nT    = length(s.eqm_dates)
    nA    = size(s.a_prob, 1)

    # Pre-compute host-age CDF for each deme (s.a_prob[:, d] is constant).
    age_cdf = cumsum(s.a_prob; dims = 1)
    for d in 1:s.nD; age_cdf[:, d] ./= age_cdf[end, d]; end

    # Pre-compute immigration CDFs (s.dProp and s.mP0 are constant).
    dProp_cdf = cumsum(s.dProp); dProp_cdf ./= dProp_cdf[end]
    mP0_cdf   = cumsum(s.mP0; dims = 1)
    for d in 1:s.nD; mP0_cdf[:, d] ./= mP0_cdf[end, d]; end

    return SimBuffer(
        zeros(Int,     s.nG, s.nD),
        zeros(Float64, s.nD, s.nL),
        zeros(Float64, s.nD, s.nL),
        zeros(Float64, s.nD, s.nL),
        zeros(Float64, s.nG, s.nD),
        zeros(Int,     s.nG * s.nD, nT),
        zeros(Float64, s.nG * s.nD, nT),
        Float64.(s.mIg0),
        Vector{Float64}(undef, s.nL),
        age_cdf,
        dProp_cdf,
        mP0_cdf,
        Vector{Int}(undef, max_n),
        Vector{Int}(undef, max_n),
        Vector{Int}(undef, max_n),
        Vector{Int}(undef, max_n),
        Vector{Float64}(undef, max_n),
        Vector{Int}(undef, max_n),
        Vector{Float64}(undef, max_n),
        Vector{Float64}(undef, max_n),
        Vector{Int}(undef, max_n),
        [Vector{Int}(undef, max_n) for _ in 1:s.nD],
        zeros(Int, s.nD),
        Vector{Int}(undef, max_n),
    )
end

"""
In-place version of simulate(): writes output to buf.result, returns true on
success or false if the population went extinct or exploded.
All major allocations are redirected into buf; the only remaining per-call
allocations are mW (nD×nD = tiny), nPer (nD = tiny), and randperm (migration).
"""
function simulate!(θ::AbstractVector{<:Real}, s::ModelState, buf::SimBuffer)
    nD = s.nD;  nV = s.nV;  nL = s.nL;  nG = s.nG

    # ── Rescale parameters ────────────────────────────────────────────────
    prop_strong = θ[1]
    f_sel       = θ[2] * 0.22
    w_sel       = θ[3] * 0.01
    migration   = θ[4] * 0.2
    co_inf      = θ[5]

    nMW = nD * (nD - 1)
    mW  = zeros(Float64, nD, nD)   # (nD×nD) — tiny; keep allocating
    pair = 1
    for d1 in 1:nD, d2 in 1:nD
        d1 == d2 && continue
        mW[d1, d2] = θ[5 + pair] * 0.1
        pair += 1
    end
    v_sel = [θ[5 + nMW + v] * 0.5 for v in 1:nV]   # nV elements — tiny
    lVs   = log1p.(-v_sel)

    # ── Selection strengths into buf.sEl ──────────────────────────────────
    fill!(buf.sEl, w_sel)
    n_strong = max(1, round(Int, prop_strong * nL))
    buf.sEl[s.sel_rank[1:n_strong]] .= f_sel

    # ── Initialise population directly into buf.gen / buf.dem ────────────
    n = 0
    for d in 1:nD
        isempty(s.tag_pre[d]) && continue
        kd  = s.kD[d]
        src = s.tag_pre[d]
        nsrc = length(src)
        for j in 1:kd
            buf.gen[n + j] = src[rand(1:nsrc)]   # sample with replacement
        end
        fill!(view(buf.dem, n+1:n+kd), d)
        n += kd
    end

    fill!(buf.result, 0)
    t_ptr = 1

    # ═══════════════════════════════════════════════════════════════════════
    for t in 1:s.n_gen

        n == 0 && return false   # population extinct

        gen_v = view(buf.gen, 1:n)
        dem_v = view(buf.dem, 1:n)

        # ── A: Genotype counts per deme (into buf.cnt) ────────────────────
        fill!(buf.cnt, 0)
        for i in 1:n;  buf.cnt[gen_v[i], dem_v[i]] += 1;  end
        nPer = vec(sum(buf.cnt; dims = 1))   # nD floats — tiny
        any(==(0), nPer) && return false

        # ★ B: NFDS fitness — all in-place, zero allocation ───────────────
        # buf.freq  = cnt' × G0 ./ nPer   →   (nD × nL)
        mul!(buf.freq, buf.cnt', s.G0)
        buf.freq ./= reshape(nPer, :, 1)
        # buf.dev = eqm_pre' .- buf.freq   →   (nD × nL)
        buf.dev .= s.eqm_pre' .- buf.freq
        # buf.adj = buf.dev .* buf.sEl'   →   (nD × nL)
        buf.adj .= buf.dev .* buf.sEl'
        # buf.nfds_fit = exp.(G0 × buf.adj')   →   (nG × nD)
        mul!(buf.nfds_fit, s.G0, buf.adj')
        buf.nfds_fit .= exp.(buf.nfds_fit)

        # Fitness per bacterium (loop — no comprehension allocation)
        fit_v = view(buf.fit, 1:n)
        for i in 1:n
            fit_v[i] = buf.nfds_fit[gen_v[i], dem_v[i]]
        end

        # ── C: Host ages — build per-deme index lists in one O(n) pass ───
        fill!(buf.n_per_d, 0)
        for i in 1:n
            d = dem_v[i]
            buf.n_per_d[d] += 1
            buf.idx_d[d][buf.n_per_d[d]] = i
        end
        h_age_v = view(buf.h_age, 1:n)
        for d in 1:nD
            nd = buf.n_per_d[d]
            nd == 0 && continue
            cdf_d = view(buf.age_cdf, :, d)
            for k in 1:nd
                i         = buf.idx_d[d][k]
                ai        = searchsortedfirst(cdf_d, rand())
                h_age_v[i] = s.a_lev[ai]
            end
        end

        # ── D: Vaccine fitness (loop — no exp.() broadcast allocation) ───
        lvp_v = view(buf.lvp, 1:n)
        fill!(lvp_v, 0.0)
        for v in 1:nV, i in 1:n
            d  = dem_v[i]
            sv = s.vstart[v, d]
            sv == Inf  && continue
            t < sv     && continue
            h_age_v[i] - (t - sv) > s.a_coh[v, d] && continue
            lvp_v[i] += s.VTmat[gen_v[i], v] * lVs[v] * s.u_pt[v, d]
        end
        for i in 1:n; fit_v[i] *= exp(lvp_v[i]); end   # in-place, no alloc

        # ── E: Offspring sampling (loops — no Poisson.() broadcast) ──────
        λ_v    = view(buf.λ,         1:n)
        offs_v = view(buf.offspring, 1:n)
        kD_sum = sum(s.kD)
        for i in 1:n
            λ_v[i] = fit_v[i] * s.kD[dem_v[i]] / nPer[dem_v[i]]
        end
        if co_inf < 1e-5
            for i in 1:n; offs_v[i] = rand(Poisson(λ_v[i])); end
        else
            r = 1.0 / co_inf
            for i in 1:n
                p = r / (r + λ_v[i])
                offs_v[i] = rand(NegativeBinomial(r, p))
            end
        end

        total_off = sum(offs_v)
        total_off == 0            && return false   # extinction
        total_off > 20 * kD_sum   && return false   # biological explosion check
        total_off > length(buf.new_gen) && return false   # buffer capacity (5× kD_sum)

        # Expand offspring into buf.new_gen / buf.new_dem
        ptr = 1
        for i in 1:n
            offs_v[i] == 0 && continue
            k = offs_v[i]
            buf.new_gen[ptr:ptr+k-1] .= gen_v[i]
            buf.new_dem[ptr:ptr+k-1] .= dem_v[i]
            ptr += k
        end
        # Zero-copy pointer swap — no data moved, only references exchanged
        buf.gen, buf.new_gen = buf.new_gen, buf.gen
        buf.dem, buf.new_dem = buf.new_dem, buf.dem
        n = total_off
        gen_v = view(buf.gen, 1:n)
        dem_v = view(buf.dem, 1:n)

        # ── F: Within-deme migration ──────────────────────────────────────
        if any(>(0), mW) && nD > 1
            perm_v = view(buf.perm, 1:n)
            perm_v .= 1:n                 # fill sequential — no allocation
            Random.shuffle!(perm_v)       # in-place Fisher-Yates — no allocation
            for i in perm_v               # was: randperm(n) → 6.4 MB/month
                d_from = dem_v[i]
                r_tot  = sum(view(mW, d_from, :))
                r_tot == 0 && continue
                rand() >= r_tot && continue
                dem_v[i] = wt_sample(view(mW, d_from, :) ./ r_tot, 1)[1]
            end
        end

        # ── G: Immigration — pre-computed CDFs, zero wt_sample allocation ─
        if migration > 0
            n_mig = rand(Poisson(migration * n))
            # Guard: immigration runs AFTER the offspring explosion check, so n can
            # grow beyond buf.gen's capacity without a separate check here.
            # Treat buffer overflow as an explosion — these parameter values produce
            # unrealistically large populations and would never be accepted as ABC particles.
            n + n_mig > length(buf.gen) && return false
            for _ in 1:n_mig
                d_to  = searchsortedfirst(buf.dProp_cdf, rand())
                g_mig = searchsortedfirst(view(buf.mP0_cdf, :, d_to), rand())
                n += 1
                buf.gen[n] = g_mig
                buf.dem[n] = d_to
            end
            gen_v = view(buf.gen, 1:n)
            dem_v = view(buf.dem, 1:n)
        end

        # ── H: Record at sampling time-points ─────────────────────────────
        if t_ptr <= length(s.eqm_dates) && t == s.eqm_dates[t_ptr]
            for i in 1:n
                buf.result[(dem_v[i]-1)*nG + gen_v[i], t_ptr] += 1
            end
            t_ptr += 1
        end

    end   # ════════════════ end generation loop ════════════════════════════

    return true
end

"""
In-place ABC distance: runs simulate! into buf, converts result to Float64
in buf.result_f, and computes JSD against the pre-converted buf.mIg0_f.
Returns the worst-case distance on simulation failure.
"""
function abc_distance!(θ::AbstractVector{<:Real}, s::ModelState, buf::SimBuffer)
    ok = simulate!(θ, s, buf)
    !ok && return Float64(length(s.eqm_dates)) * log(2)
    buf.result_f .= buf.result
    return total_jsd(buf.result_f, buf.mIg0_f)
end
