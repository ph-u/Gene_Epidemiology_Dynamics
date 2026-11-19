# preprocess--mass.jl — convert raw mass.input to the v2 data format
# Port of R's preprocess--mass.r.
#
# Input:   data/mass.input      (raw Massachusetts genomics data, 616 × 5430)
# Output:  data/pp-mass.tsv     (v2 format: Time|Deme|Age|VT1|SC|genes…)
#          data/pp-roll.tsv     (vaccine rollout schedule)
#          data/input.csv       (updated file paths)
#
# Run from the repo root:
#   julia src/preprocess--mass.jl

using CSV, DataFrames, Random, Statistics

Random.seed!(20261029)

# ─── Read inputs ───────────────────────────────────────────────────────────────
cfg = CSV.read("data/input.csv", DataFrame)
raw = CSV.read("data/mass.input", DataFrame; delim = '\t')
println("Read $(nrow(raw)) isolates, $(ncol(raw)) columns")
println("Time distribution: ", sort([(t, count(==(t), raw.Time)) for t in unique(raw.Time)]))
println("VT  distribution: ", sort([(v, count(==(v), raw.VT))   for v in unique(raw.VT)]))
println("SC unique values:  $(length(unique(raw.SC)))")

# ─── Build v2 data table ───────────────────────────────────────────────────────
# Required column order: Time | Deme | Age | VT1 | SC | <gene columns>
# "SC" must be the LAST metadata column before the gene columns.
sc_pos    = findfirst(==("SC"), names(raw))
gene_cols = names(raw)[sc_pos+1:end]   # pbp*, CLS*, UNCLS

# Age imputation:
# The original dataset (Huang et al. 2009, PIDJ) surveyed children aged 3–84 months.
# Age is not in the public data; we impute uniformly over that range and then
# shift by Time (months since vaccine introduction), capping at 84 months.
# This is identical to the R implementation: set.seed(20261029); sample(3:84) + Time.
raw_ages = rand(3:84, nrow(raw))             # uniform draw from 3..84
raw.Age  = clamp.(raw_ages .+ raw.Time, 3, 84)

# Assemble v2 columns: keep all gene columns in their original order.
pp = hcat(
    DataFrame(Time = raw.Time,
              Deme = "d1",        # single location (Massachusetts)
              Age  = raw.Age,
              VT1  = raw.VT,      # renamed VT → VT1  (v2 naming convention)
              SC   = raw.SC),
    raw[:, gene_cols]
)
CSV.write("data/pp-mass.tsv", pp; delim = '\t')
println("\nWritten $(nrow(pp)) rows × $(ncol(pp)) cols → data/pp-mass.tsv")

# ─── Vaccine rollout schedule ──────────────────────────────────────────────────
# PCV7 introduced in Massachusetts at Month 0.
# Catch-up programme covered children ≤59 months (Lexau et al. 2005).
# AgeWindow = maximum observed age covers all children in the study cohort,
# matching the R preprocess behaviour.
roll = DataFrame(
    Deme       = ["d1"],
    Vaccine    = ["VT1"],
    StartMonth = [0],
    AgeWindow  = [maximum(pp.Age)],
    Uptake     = [1.0]
)
CSV.write("data/pp-roll.tsv", roll; delim = '\t')
println("Written rollout → data/pp-roll.tsv")

# ─── Update input.csv with new file paths ─────────────────────────────────────
cfg[cfg.Type .== "data",    :Value] .= ["../data/pp-mass.tsv"]
cfg[cfg.Type .== "rollout", :Value] .= ["../data/pp-roll.tsv"]
CSV.write("data/input.csv", cfg)
println("Updated data/input.csv")
println("\nDone.  Run model with:  julia --threads auto src/model.jl ../data/input.csv 1")
