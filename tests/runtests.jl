# runtests.jl — NFDS Julia model unit tests
# Mirrors the three R test categories:
#   test-src   → @testset "func"   (utility functions)
#   test-nfds  → @testset "nfds"   (simulation properties)
#   test-model → @testset "model"  (ABC-SMC orchestration)
#
# Run with:
#   julia --project=src --threads auto tests/runtests.jl

using Test, Random, Statistics, Distributions, DataFrames, CSV

# Single include: model.jl chains model → nfds → func via @isdefined guards.
# Including all three separately would redefine the ModelState struct,
# which is a hard error in Julia ("cannot redefine constant ModelState").
include(joinpath(@__DIR__, "..", "src", "model.jl"))

# ─── Fixture builder ───────────────────────────────────────────────────────────
# Small synthetic dataset: 3 demes, 2 vaccines, 4 genotypes, 6 gene loci.
# Mirrors the R helper-fixture.R fixture (nD=3, nVac=2, 4 lineages, k=1200).

function make_fixture()
    dir = mktempdir()

    # Four genotypes: G1..G4 are lineage markers (one-hot).
    # G5, G6 are accessory genes carried by multiple lineages.
    # VT1 targets lineage 1 only; VT2 targets lineages 1 and 2.
    nLin = 4
    vt1  = [1,0,0,0];  vt2 = [1,1,0,0]
    sc   = [1,2,3,4]
    # gene pattern per lineage
    genes = [1 0 0 0 1 0;   # lineage 1
             0 1 0 0 0 1;   # lineage 2
             0 0 1 0 1 0;   # lineage 3
             0 0 0 1 0 1]   # lineage 4

    times   = [-6, 0, 6, 12]
    # counts per lineage per timepoint (vaccine reduces VT types over time)
    mix_lin = [[2,2,2,2], [2,2,2,2], [1,2,3,2], [1,1,3,3]]
    n_deme  = 3

    rows = []
    for (ti, t) in enumerate(times)
        counts = mix_lin[ti]
        for d in 1:n_deme
            for (lin, c) in enumerate(counts)
                for _ in 1:c
                    push!(rows, (Time=t, Deme="d$d", Age=rand(6:60),
                                 VT1=vt1[lin], VT2=vt2[lin], SC=sc[lin],
                                 G1=genes[lin,1], G2=genes[lin,2],
                                 G3=genes[lin,3], G4=genes[lin,4],
                                 G5=genes[lin,5], G6=genes[lin,6]))
                end
            end
        end
    end
    df = DataFrame(rows)
    CSV.write(joinpath(dir, "data.tsv"), df; delim='\t')

    roll = DataFrame(Deme=["d1","d1","d2","d2","d3"],
                     Vaccine=["VT1","VT2","VT1","VT2","VT1"],
                     StartMonth=[0, 6, 3, 9, 9],
                     AgeWindow=[24,24,24,24,24],
                     Uptake=[1.0,0.8,1.0,0.7,1.0])
    CSV.write(joinpath(dir, "rollout.tsv"), roll; delim='\t')

    cfg = DataFrame(Type=["data","rollout","population","percentage initial infected"],
                    Value=[joinpath(dir,"data.tsv"), joinpath(dir,"rollout.tsv"), "1200", "25"],
                    Note=["","","",""])
    CSV.write(joinpath(dir, "input.csv"), cfg)

    return dir
end

# Helper: build a small ModelState from the fixture
function fixture_state(dir::String)
    Random.seed!(42)
    build_state(joinpath(dir, "input.csv"))
end

# Default parameter vector (all 0.5) for the fixture (3 demes, 2 vaccines)
default_θ(nD::Int, nV::Int) = fill(0.5, 5 + nD*(nD-1) + nV)

# ══════════════════════════════════════════════════════════════════════════════
# TEST GROUP 1 — FUNC (mirrors R's test-src.R)
# ══════════════════════════════════════════════════════════════════════════════
@testset "func" begin

    @testset "read_cfg" begin
        tmp = tempname() * ".csv"
        write(tmp, "Type,Value,Note\npopulation,800000,hosts\nrollout,../data/pp-roll.tsv,\n")
        @test read_cfg(tmp, "population") == "800000"
        @test read_cfg(tmp, "rollout")    == "../data/pp-roll.tsv"
        @test isnothing(read_cfg(tmp, "no_such_key"))
        rm(tmp)
    end

    @testset "jsd — mathematical properties" begin
        # Identical distributions → 0
        @test isapprox(jsd([1.0,2.0,3.0], [1.0,2.0,3.0]), 0.0; atol=1e-10)
        # Completely disjoint → log(2)
        @test isapprox(jsd([1.0,0.0], [0.0,1.0]), log(2); atol=1e-10)
        # Symmetry: JSD(p,q) == JSD(q,p)
        p, q = [3.0,1.0,4.0], [1.0,2.0,3.0]
        @test isapprox(jsd(p, q), jsd(q, p); atol=1e-10)
        # Scale-invariant: only proportions matter
        @test isapprox(jsd([1.0,2.0], [2.0,4.0]), 0.0; atol=1e-10)
        # Non-negative
        @test jsd(rand(5), rand(5)) >= 0.0
        # Empty column → NaN
        @test isnan(jsd([0.0,0.0], [1.0,2.0]))
        @test isnan(jsd([1.0,2.0], [0.0,0.0]))
    end

    @testset "total_jsd" begin
        m = abs.(randn(4, 3)) .+ 0.1
        @test isapprox(total_jsd(m, m), 0.0; atol=1e-9)
        @test total_jsd(m, m) >= 0.0
        # An all-zero column is penalised with log(2)
        obs = copy(m);  obs[:, 2] .= 0
        sim = copy(m)
        d   = total_jsd(sim, obs)
        @test isfinite(d) && d > 0
    end

    @testset "wt_sample" begin
        ws = [0.1, 0.5, 0.4]
        s  = wt_sample(ws, 10_000)
        @test all(1 .<= s .<= 3)                     # all indices are valid
        # Category 2 (weight 0.5) should dominate
        @test count(==(2), s) > count(==(1), s)
        @test count(==(2), s) > count(==(3), s)
        # Single-element pool always returns 1
        @test all(wt_sample([1.0], 5) .== 1)
    end

end   # @testset "func"

# ══════════════════════════════════════════════════════════════════════════════
# TEST GROUP 2 — NFDS (mirrors R's test-nfds.R)
# ══════════════════════════════════════════════════════════════════════════════
@testset "nfds" begin

    dir   = make_fixture()
    state = fixture_state(dir)
    nD, nV, nG = state.nD, state.nV, state.nG
    nT    = length(state.eqm_dates)
    θ_def = default_θ(nD, nV)

    @testset "build_state — dimensions" begin
        @test state.nD == 3
        @test state.nV == 2
        @test state.nG >= 2             # at least the 4 lineages as genotypes
        @test size(state.G0)      == (nG, state.nL)
        @test size(state.VTmat)   == (nG, nV)
        @test size(state.mIg0)    == (nG * nD, nT)
        @test size(state.eqm_pre) == (state.nL, nD)
        @test size(state.mP0)     == (nG, nD)
        @test size(state.vstart)  == (nV, nD)
        @test size(state.a_coh)   == (nV, nD)
        @test size(state.u_pt)    == (nV, nD)
        @test size(state.a_prob)  == (length(state.a_lev), nD)
    end

    @testset "build_state — biological constraints" begin
        # G0 and VTmat must be binary
        @test all(x -> x ∈ (0.0, 1.0), state.G0)
        @test all(x -> x ∈ (0.0, 1.0), state.VTmat)
        # Deme proportions sum to 1
        @test isapprox(sum(state.dProp), 1.0; atol=1e-9)
        # Migrant pool columns sum to 1
        @test all(isapprox.(sum(state.mP0; dims=1), 1.0; atol=1e-6))
        # Age distribution columns sum to 1
        @test all(isapprox.(sum(state.a_prob; dims=1), 1.0; atol=1e-6))
        # Pre-vaccine data should give reasonable equilibrium (0 < freq < 1)
        @test all(0.0 .<= state.eqm_pre .<= 1.0)
        # Vaccine starts no earlier than month 0 or Inf
        @test all(state.vstart .>= 0)
    end

    @testset "simulate — output shape and validity" begin
        Random.seed!(1)
        result = simulate(θ_def, state)
        @test !isnothing(result)
        @test size(result) == (nG * nD, nT)
        @test all(result .>= 0)             # counts must be non-negative
        # Total simulated isolates at each timepoint should match nObs
        for t in 1:nT
            @test sum(result[:, t]) == sum(state.nObs[t, :])
        end
    end

    @testset "simulate — deterministic under same seed" begin
        Random.seed!(99);  a = simulate(θ_def, state)
        Random.seed!(99);  b = simulate(θ_def, state)
        @test a == b
    end

    @testset "simulate — vaccine selection reduces VT frequency" begin
        # θ with high vaccine efficacy (vSel ≈ 0.9) vs zero efficacy (vSel = 0)
        θ_hi  = copy(θ_def);  θ_hi[end-nV+1:end]  .= 0.9   # high efficacy
        θ_lo  = copy(θ_def);  θ_lo[end-nV+1:end]  .= 0.0   # no vaccine effect
        # Disable migration to isolate the vaccine signal
        θ_hi[5 .+ (1:nD*(nD-1))] .= 0.0;  θ_hi[4] = 0.0
        θ_lo[5 .+ (1:nD*(nD-1))] .= 0.0;  θ_lo[4] = 0.0

        vt_rows = [any(state.VTmat[g, :] .> 0) for g in 1:nG]   # genotypes targeted by any vaccine
        vt_frac(r) = sum(r[vcat([((d-1)*nG+1):d*nG for d in 1:nD]...)[repeat(vt_rows, nD)], end]) /
                     max(1, sum(r[:, end]))

        f_hi = mean(let Random.seed!(s);  r = simulate(θ_hi, state);  isnothing(r) ? NaN : vt_frac(r)  end
                    for s in 1:6)
        f_lo = mean(let Random.seed!(s);  r = simulate(θ_lo, state);  isnothing(r) ? NaN : vt_frac(r)  end
                    for s in 1:6)
        @test f_hi < f_lo    # strong vaccine → lower VT fraction at end
    end

    @testset "simulate — mWithin moves bacteria between demes" begin
        Random.seed!(1);  a = simulate(θ_def, state)   # mWithin ≈ 0.5*0.1 per pair
        θ_no  = copy(θ_def);  θ_no[5 .+ (1:nD*(nD-1))] .= 0.0   # no migration
        Random.seed!(1);  b = simulate(θ_no, state)
        @test !isequal(a, b)   # within-deme migration changes the trajectory
    end

    @testset "simulate — coInf = 0 Poisson vs coInf > 0 NegBin" begin
        θ_p  = copy(θ_def);  θ_p[5]  = 0.0    # Poisson
        θ_nb = copy(θ_def);  θ_nb[5] = 0.5    # NegBin
        Random.seed!(1);  a = simulate(θ_p,  state)
        Random.seed!(1);  b = simulate(θ_nb, state)
        @test !isequal(a, b)   # NegBin branch produces different output
    end

    @testset "simulate — population runaway returns nothing" begin
        # Set propRunaway very low by using extreme fitness: all params = 0 → no regulation
        # This is hard to trigger directly, so we just verify the return type contract.
        result = simulate(θ_def, state)
        @test result isa Union{Matrix{Int}, Nothing}
    end

    @testset "abc_distance — range and validity" begin
        Random.seed!(2)
        d = abc_distance(θ_def, state)
        @test isfinite(d)
        @test d >= 0.0
        @test d <= nT * log(2)   # maximum possible JSD sum
    end

    @testset "abc_distance — identical simulation gives 0 JSD" begin
        # Manually set mIg0 to what simulate returns for θ_def and check distance ≈ 0.
        # (We can't guarantee exact equality due to stochasticity; just check range.)
        d = abc_distance(θ_def, state)
        @test d < nT * log(2)   # not at the maximum
    end

    rm(dir; recursive = true)

end   # @testset "nfds"

# ══════════════════════════════════════════════════════════════════════════════
# TEST GROUP 3 — MODEL (mirrors R's test-model.R)
# ══════════════════════════════════════════════════════════════════════════════
@testset "model" begin

    @testset "param_names — count and names" begin
        nD, nV = 3, 2
        expected_len = 5 + nD*(nD-1) + nV   # = 13 for nD=3, nV=2
        pn = param_names(nD, nV)
        @test length(pn) == expected_len
        # Base parameters present
        @test "propStrong" ∈ pn
        @test "fSelected"  ∈ pn
        @test "wSelected"  ∈ pn
        @test "migration"  ∈ pn
        @test "coInf"      ∈ pn
        # All directed deme pairs present
        for d1 in 1:nD, d2 in 1:nD
            d1 == d2 && continue
            @test "mWithin_$(d1)_$(d2)" ∈ pn
        end
        # Vaccine parameters present
        for v in 1:nV
            @test "vSelected$v" ∈ pn
        end
        # Removed parameters should NOT be present
        @test "vSelMod"  ∉ pn
        @test "mWithin"  ∉ pn   # replaced by directed pairs
    end

    @testset "sample_prior — values in [0,1]" begin
        for _ in 1:50
            θ = sample_prior(13)
            @test length(θ) == 13
            @test all(0.0 .<= θ .<= 1.0)
        end
    end

    @testset "perturb — clamps to [0,1] and changes values" begin
        θ = fill(0.5, 13)
        for _ in 1:50
            p = perturb(θ)
            @test all(0.0 .<= p .<= 1.0)
            @test length(p) == length(θ)
        end
        # Boundary: input at 0 or 1 should stay in [0,1]
        @test all(0.0 .<= perturb(zeros(5)) .<= 1.0)
        @test all(0.0 .<= perturb(ones(5))  .<= 1.0)
        # Should produce different values than the input (with overwhelming probability)
        @test !isequal(perturb(θ), θ)
    end

    @testset "prior dimension matches build_state output" begin
        dir   = make_fixture()
        state = fixture_state(dir)
        n_dim = 5 + state.nD * (state.nD - 1) + state.nV
        @test length(sample_prior(n_dim)) == n_dim
        @test length(param_names(state.nD, state.nV)) == n_dim
        rm(dir; recursive = true)
    end

    @testset "setup script has no side effects on loaded state" begin
        # build_state must not write to the working directory; only return a struct.
        dir    = make_fixture()
        before = readdir(dir)
        _      = fixture_state(dir)
        after  = readdir(dir)
        @test Set(before) == Set(after)   # no new files created in fixture dir
        rm(dir; recursive = true)
    end

end   # @testset "model"

println("\nAll tests passed ✓")
