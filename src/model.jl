# model.jl — ABC-SMC parameter estimation
# Equivalent to R's model.r + BRREWABC, but all parallelism runs in shared memory.
#
# WHY JULIA THREADS INSTEAD OF BRREWABC:
#   BRREWABC communicates between the master process and worker processes by
#   writing temporary CSV files to disk.  Profiling showed this consumed 98.8%
#   of wall time.  Threads.@threads evaluates proposals in parallel within ONE
#   process — no files, no serialisation, no process-launch overhead.
#   Computational effect: N threads → N simultaneous simulations, same as
#   N separate R processes but with 100× less per-particle overhead.
#
# Usage (from the repo root):
#   julia --project=src --threads auto src/model.jl ../data/input.csv 1
#
#   --threads auto    = use all available physical CPU cores
#   second argument   = integer seed (matches the R seed convention)

# Same pattern as nfds.jl: guard against double-include and use @__DIR__ for
# an absolute path that works from any working directory.
@isdefined(abc_distance) || include(joinpath(@__DIR__, "nfds.jl"))

# CSV and DataFrames are already loaded by nfds.jl, but listing them here
# makes the dependency explicit and avoids confusion for readers of model.jl.
using CSV, DataFrames, Random, Statistics
using Printf: @printf

# ─── Parameter layout (matches R model.r prior) ───────────────────────────────
# θ is a vector of length  5 + nD*(nD-1) + nV  with all values in [0, 1].
# Biological rescaling to real units happens inside simulate() (see nfds.jl).
#
# Idx  Name              Biological range (after rescaling inside simulate)
# ───  ────────────────  ────────────────────────────────────────────────────
#  1   propStrong        fraction of loci under strong NFDS   (0 → 1)
#  2   fSelected         strong-locus NFDS strength            (0 → 0.22 / month)
#  3   wSelected         weak-locus  NFDS strength             (0 → 0.01 / month)
#  4   migration         external immigration rate             (0 → 0.2  / month)
#  5   coInf             offspring overdispersion (0 = Poisson)
#  6…  mWithin_d1_d2     one per directed deme pair            (0 → 0.1  / month)
#  …   vSelected1…V      per-vaccine efficacy                  (0 → 0.5  / month)

"""Return the ordered parameter names for a model with nD demes and nV vaccines."""
function param_names(nD::Int, nV::Int)
    base = ["propStrong", "fSelected", "wSelected", "migration", "coInf"]
    mw   = ["mWithin_$(d1)_$(d2)" for d1 in 1:nD for d2 in 1:nD if d1 ≠ d2]
    vs   = ["vSelected$v" for v in 1:nV]
    return [base; mw; vs]
end

"""Sample one parameter vector uniformly from [0,1]^d (the flat prior)."""
sample_prior(d::Int) = rand(Float64, d)

"""Perturb a particle with Gaussian noise, then clamp strictly to [0,1].
σ = 0.1 matches the typical BRREWABC Gaussian kernel width."""
perturb(θ::AbstractVector, σ::Float64 = 0.1) =
    clamp.(θ .+ σ .* randn(length(θ)), 0.0, 1.0)

# ─── ABC-SMC loop ─────────────────────────────────────────────────────────────

"""
Run ABC-SMC.  Parameters match the R+BRREWABC implementation:
  pilot_n     = 200   — pilot simulations to establish the initial threshold
  n_particles = 500   — accepted particles per SMC generation
  n_smc_gen   = 20    — number of successive tightening rounds

Checkpoints are written to `../data/run_<seed>/smc_gen***.csv` after each
generation.  The final posterior is also saved as `posterior.csv`.
"""
function run_abc_smc(state::ModelState, seed::Int;
                     pilot_n     ::Int = 200,
                     n_particles ::Int = 500,
                     n_smc_gen   ::Int = 20)

    Random.seed!(seed)
    n_dim   = 5 + state.nD * (state.nD - 1) + state.nV
    pnames  = param_names(state.nD, state.nV)
    out_dir = mkpath("../data/run_$seed")
    n_cpu   = Threads.nthreads()
    println("ABC-SMC | seed=$seed | $n_cpu threads | $n_dim parameters\n")

    # ── Phase 1: pilot run ────────────────────────────────────────────────
    # Draw pilot_n parameter vectors from the flat prior, evaluate in parallel,
    # and use the median distance as the starting acceptance threshold.
    # This matches the BRREWABC PILOT_N step.
    println("Phase 1 — pilot (n=$pilot_n, using $n_cpu threads)…")
    pilot_θ = [sample_prior(n_dim) for _ in 1:pilot_n]
    pilot_d = Vector{Float64}(undef, pilot_n)
    Threads.@threads for i in 1:pilot_n     # ← N proposals evaluated simultaneously
        pilot_d[i] = abc_distance(pilot_θ[i], state)
    end

    finite_d = filter(isfinite, pilot_d)
    isempty(finite_d) && error("All pilot simulations failed — check data/input.csv paths")
    threshold = quantile(finite_d, 0.5)   # accept best 50% of pilot runs

    ok        = findall(d -> isfinite(d) && d <= threshold, pilot_d)
    particles = [pilot_θ[i] for i in ok]
    dists     = [pilot_d[i]  for i in ok]
    # Top up to n_particles if acceptance rate was low
    while length(particles) < n_particles
        θ = sample_prior(n_dim)
        d = abc_distance(θ, state)
        d <= threshold && (push!(particles, θ); push!(dists, d))
    end
    particles = particles[1:n_particles];  dists = dists[1:n_particles]
    @printf("  ε₀ = %.4f  (accepted %d / %d pilot runs)\n\n", threshold, length(ok), pilot_n)

    # ── Phase 2: successive SMC generations ───────────────────────────────
    # Each generation:
    #   (a) Set ε = median of current accepted distances (halves roughly each round).
    #   (b) Generate n_particles × 20 perturbed proposals.
    #   (c) Evaluate all in parallel with Threads.@threads — zero file I/O.
    #   (d) Accept the best n_particles.
    #   (e) Save CSV checkpoint.
    for g in 1:n_smc_gen
        threshold = quantile(dists, 0.5)

        n_cand = n_particles * 20                      # oversample heavily
        cands  = [perturb(rand(particles)) for _ in 1:n_cand]
        cand_d = Vector{Float64}(undef, n_cand)
        Threads.@threads for i in 1:n_cand             # ← the parallelism
            cand_d[i] = abc_distance(cands[i], state)
        end

        ok = findall(d -> isfinite(d) && d <= threshold, cand_d)
        if length(ok) < n_particles
            @printf("  Gen %2d | ε=%.4f | only %d/%d accepted — no update\n",
                    g, threshold, length(ok), n_cand)
            continue
        end
        top       = ok[sortperm(cand_d[ok])[1:n_particles]]
        particles = [cands[i] for i in top]
        dists     = cand_d[top]

        # ── Checkpoint: one CSV row per particle ─────────────────────────
        mat = reduce(hcat, particles)'       # (n_particles × n_dim)
        df  = DataFrame(mat, pnames)
        df.distance = dists
        path = "$out_dir/smc_gen$(lpad(g, 3, '0')).csv"
        CSV.write(path, df)
        @printf("  Gen %2d | ε=%.4f | accept=%d/%d | mean_d=%.4f | %s\n",
                g, threshold, length(ok), n_cand, mean(dists), basename(path))
    end

    # ── Final posterior ───────────────────────────────────────────────────
    final_df = DataFrame(reduce(hcat, particles)', pnames)
    final_df.distance = dists
    CSV.write("$out_dir/posterior.csv", final_df)
    println("\nFinal posterior  → $out_dir/posterior.csv")
    println("Checkpoints      → $out_dir/smc_gen***.csv")
end

# ─── Entry point ───────────────────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) < 2 && error("""
    Usage:
      julia --project=src --threads auto src/model.jl <input_csv> <seed>

    Example (10 cores, seed 1):
      julia --project=src --threads 10 src/model.jl ../data/input.csv 1
    """)

    input_csv  = ARGS[1]
    seed       = parse(Int, ARGS[2])
    state_path = "../data/run_$seed/setup_state.jld2"

    if isfile(state_path)
        println("Loading cached state from $state_path…")
        state = load_state(state_path)
    else
        println("Building model state from $input_csv…")
        mkpath(dirname(state_path))
        state = build_state(input_csv)
        save_state(state_path, state)
        println("State saved → $state_path")
    end

    run_abc_smc(state, seed)
end
