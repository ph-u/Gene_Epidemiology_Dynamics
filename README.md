# NFDS Population Genetics Model — v3 Documentation

> **Model**: Negative-Frequency-Dependent Selection (NFDS) + vaccine pressure  
> **Method**: Approximate Bayesian Computation – Sequential Monte Carlo (ABC-SMC)  
> **Language**: Julia 1.11  
> **Reference**: Corander et al. (2017) *Nature Ecology & Evolution* doi:10.1038/s41559-017-0337-x

---

## 1. Prerequisites

Install Julia 1.11 from https://julialang.org/downloads/ and verify:

```bash
julia --version   # should print 1.11.x
```

Install all Julia packages (once, from the repo root):

```bash
julia --project=src -e 'using Pkg; Pkg.instantiate()'
```

---

## 2. Data preparation

The raw `data/mass.input` must be converted to the v2 format before fitting:

```bash
# From the repo root:
julia src/preprocess--mass.jl
```

**What this does:**

| step | detail |
|---|---|
| Adds `Deme` column | constant `"d1"` (single location) |
| Adds `Age` column | imputed 3–84 months (Huang et al. 2009); shift by Time, cap at 84 |
| Renames `VT` → `VT1` | v2 naming convention for multi-vaccine support |
| Writes `data/pp-mass.tsv` | tab-separated, v2 column order |
| Writes `data/pp-roll.tsv` | PCV7 rollout schedule (start month 0, all children) |
| Updates `data/input.csv` | points to new file paths |

**v2 column order required by `build_state()`:**

```
Time | Deme | Age | VT1 [| VT2 …] | SC | <gene columns>
```

`SC` must be the **last** metadata column immediately before the gene columns.

---

## 3. Running the model

```bash
# Standard run on the HPC (N = number of CPU cores available):
julia --project=src --threads N src/model.jl ../data/input.csv <seed>

# Example: 10 cores, seed 1
julia --project=src --threads 10 src/model.jl ../data/input.csv 1

# Multiple seeds in parallel (each is an independent chain):
for seed in 1 2 3 4 5; do
  julia --project=src --threads 10 src/model.jl ../data/input.csv $seed &
done
wait
```

**Why threads instead of separate processes (R+BRREWABC)?**  
The R implementation used BRREWABC which wrote every particle evaluation to a
temporary CSV file and read it back. Profiling showed this consumed **98.8% of
wall time**. Julia threads share the loaded `ModelState` in RAM. No files are
written per particle, so overhead per evaluation is negligible.

**ABC-SMC parameters (matching the R run):**

| parameter | value | meaning |
|---|---|---|
| `pilot_n` | 200 | pilot simulations to find the initial threshold ε₀ |
| `n_particles` | 500 | accepted particles per SMC generation |
| `n_smc_gen` | 20 | number of successive tightening rounds |
| `σ` | 0.1 | Gaussian perturbation kernel width |

To change these defaults, edit `run_abc_smc(...)` in `src/model.jl`.

---

## 4. Output files

All outputs land in `data/run_<seed>/`:

| file | content |
|---|---|
| `setup_state.jld2` | pre-computed model state (cached; reused on restart) |
| `smc_gen001.csv` … `smc_gen020.csv` | 500 particles × (n_params + 1) after each SMC generation |
| `posterior.csv` | final accepted particles from the last generation |

Each CSV row is one accepted particle. Columns are the parameter names followed by `distance` (the JSD value that made it accepted).

---

## 5. Interpreting the output

### 5.1 Parameter names and biological meaning

| parameter | rescaled range | interpretation |
|---|---|---|
| `propStrong` | 0 → 1 | fraction of gene loci under **strong** NFDS |
| `fSelected` | 0 → 0.22/month | NFDS strength for strong-locus genes |
| `wSelected` | 0 → 0.01/month | NFDS strength for weak-locus genes |
| `migration` | 0 → 0.2/month | rate of external immigration into the host population |
| `coInf` | ≥ 0 | co-infection overdispersion; 0 = Poisson offspring |
| `mWithin_d1_d2` | 0 → 0.1/month | migration rate from deme d1 to deme d2 |
| `vSelected1` … | 0 → 0.5/month | per-vaccine fitness reduction on targeted bacteria |

The raw `θ` values stored in the CSV are in **[0, 1]**. To convert to biological
units, apply the same rescaling as `simulate()` in `src/nfds.jl` (multiply by
the range shown above).

### 5.2 Reading posterior CSVs in Julia

```julia
using CSV, DataFrames, Statistics

posterior = CSV.read("data/run_1/posterior.csv", DataFrame)

# Point estimates (posterior medians)
params = ["propStrong","fSelected","wSelected","migration","coInf","vSelected1"]
for p in params
    println("$p: $(round(median(posterior[!, p]); digits=4))")
end

# 95% credible intervals
for p in params
    lo, hi = quantile(posterior[!, p], [0.025, 0.975])
    println("$p: [$( round(lo;digits=3)), $(round(hi;digits=3))]")
end
```

### 5.3 Convergence diagnostics

**SMC generation summaries** are printed to the console during the run. Look for:
- `mean_d` decreasing steadily across generations → the posterior is tightening
- `accept = k / 10000` ratio: should stay above ~1% throughout; if it drops near
  zero before `n_smc_gen` is reached, the model may be over-constrained

**Across seeds**: load `posterior.csv` from multiple seeds and compare marginal
distributions. Use `coda` (R) or `MCMCChains.jl` (Julia) for Gelman-Rubin R̂:

```julia
# In Julia, using MCMCChains.jl:
using MCMCChains
chains = Chains(cat([Matrix(CSV.read("data/run_$s/posterior.csv", DataFrame)[:, 1:end-1])
                     for s in 1:5]...; dims=3))
summarystats(chains)   # R̂ ≤ 1.01 indicates good mixing across seeds
```

### 5.4 Biological interpretation of key parameters

| observation | interpretation |
|---|---|
| `propStrong` close to 1 | most loci contribute to NFDS → strong diversity maintenance |
| `vSelected1` near 0.5 | maximum vaccine pressure; VT types experience ~40% fitness reduction |
| `coInf` > 0.1 | bacterial colonisation is bursty (super-spreading); NegBin offspring |
| `mWithin_d1_d2` > 0.05 | substantial bacterial flow between sub-populations |
| High `migration` | significant influx from outside the study population |

---

## 6. Quick restart after interruption

If the run was interrupted, `model.jl` automatically loads the cached
`setup_state.jld2` on the next invocation — no data re-processing needed:

```bash
julia --project=src --threads 10 src/model.jl ../data/input.csv 1
# prints: "Loading cached state from ../data/run_1/setup_state.jld2…"
```

To force a full rebuild (e.g. after changing the data):

```bash
rm data/run_1/setup_state.jld2
julia --project=src --threads 10 src/model.jl ../data/input.csv 1
```

---

## 7. Directory structure

```
├── src/
│   ├── func.jl              # utility functions (JSD, config reader, sampler)
│   ├── nfds.jl              # model state build + NFDS forward simulation
│   ├── model.jl             # ABC-SMC with Threads.@threads parallelism
│   ├── preprocess--mass.jl  # raw mass.input → v2 TSV format
│   └── Project.toml         # Julia package manifest (pinned to Julia 1.11)
├── tests/
│   └── runtests.jl          # unit tests (func / nfds / model categories)
├── data/
│   ├── input.csv            # master config (updated by preprocess--mass.jl)
│   ├── mass.input           # raw Massachusetts genomics data (read-only)
│   ├── pp-mass.tsv          # v2-format data  (created by preprocess--mass.jl)
│   ├── pp-roll.tsv          # vaccine rollout (created by preprocess--mass.jl)
│   └── run_<seed>/          # output per seed
│       ├── setup_state.jld2
│       ├── smc_gen***.csv
│       └── posterior.csv
├── raw/
│   └── seed.csv             # hash numbers for numeric reproducibility of simulations
├── Dockerfile
└── env.yaml
```

## AI Disclaimer

All code (ver 2 and 3) and documentation produced for this pipeline was initially generated by Claude, but reviewed by humans.
