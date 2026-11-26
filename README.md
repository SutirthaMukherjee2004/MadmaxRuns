# MadmaxRuns

Julia implementation of the analysis chain of *Signal-agnostic search for narrow lines on
unknown backgrounds* (SQUARES test statistic with Fourier whitening). It also contains the
studies that reproduce every figure and table of the paper, a CairoMakie plotting package and
PBS batch scripts.

## File tree

```
MadmaxRuns/
├── Project.toml
├── README.md
├── src/
│   ├── MadmaxRuns.jl            module, includes and exports
│   ├── model/
│   │   ├── grid.jl              frequency grid (25 000 bins of 2 kHz from 4.218 MHz)
│   │   ├── background.jl        background b(f), Eq. (bkg), r_k ~ N(0,1)
│   │   ├── noise.jl             σ_i = σ_0 + σ_1 b_i, homo-/heteroscedastic, Eq. (C1)
│   │   └── lineshapes.jl        Gaussian, Lorentzian, top-hat, axion (SHM) lines, Eq. (sig)
│   ├── filter/
│   │   ├── savitzky_golay.jl    SG / SGW kernels, residual kernel g, G(f), r = G_lin x
│   │   └── weighted_sg.jl       SG weighted by w_i = 1/σ̂_i² (FFT correlations + per-bin solve)
│   ├── noise/estimate.jl        windowed fit (Eq. s7), LSQ fit (Eq. s6), iterative estimate
│   ├── whitening/fourier.jl     mode set M, κ, division by G(f) on M
│   ├── boundary/
│   │   ├── boundary.jl          spike responses, M̃, SVD, removal z' = z − Σ (v_kᵀz) v_k
│   │   └── tau0.jl              τ0 exact (Eqs. f5–f7) and closed form (Eq. f10)
│   ├── stats/
│   │   ├── special.jl           Ψ_q, Ψ̄_q, Φ, Φ⁻¹, tilted χ² tail moments
│   │   ├── bivariate_normal.jl  Φ₂ (Gauss–Legendre)
│   │   ├── squares.jl           T from RunStatistics.jl, P_BC (RunStatistics.jl)
│   │   ├── pk.jl                P_K, Eq. (g5): E_q first order and second order (Eq. g7), tables
│   │   └── pk_exact.jl          exact E_q, Eq. (g6), by Monte Carlo
│   ├── pipeline/
│   │   ├── config.jl            PipelineConfig (n, m, kind, τ, W, estimator, tol, K, P_K order, …)
│   │   ├── setup.jl             FilterSetup: everything that depends only on filter and τ
│   │   └── chain.jl             run_chain (steps 2–10), simulate_spectrum (step 1)
│   ├── nulls/nulls.jl           H0^iid and H0^filt Monte-Carlo nulls, MC p values
│   ├── io/                      logging, results (.jls/.csv), thread settings, command line
│   └── studies/                 one file per figure/table group, registry.jl lists them
├── scripts/
│   ├── setup.jl                 installs FFTW, SpecialFunctions, RunStatistics (and CairoMakie)
│   ├── run_pipeline.jl          the chain for one SG filter (n, m) on one spectrum
│   ├── run_study.jl             one study (or all) from the command line
│   ├── selftest.jl              stepwise tests + quick run of every study
│   ├── merge_tables.jl          merges the table sub-jobs, writes table_shapes4.tex
│   └── latex_table.jl
├── batch/
│   ├── config.sh                paths, queue, cores, memory, walltime
│   ├── jobs_stage1.txt          thresholds τ_min, τ0 and independent studies
│   ├── jobs_stage2.txt          studies that use τ_min, τ0; table split by shape and width
│   ├── pbs_array.sh             one PBS sub-job = one line of a job list
│   ├── pbs_single.sh, pbs_finish.sh
│   ├── pbs_all.sh               submits everything with dependencies
│   └── run_local.sh             the same on one machine without a scheduler
├── plots/
│   ├── Project.toml             MadmaxRunsPlots (CairoMakie)
│   ├── src/                     theme and one function per paper figure
│   └── scripts/make_figures.jl
└── test/                        stepwise tests of every step (runtests.jl)
```

## Flowchart → code

| Step | What | Code |
|---|---|---|
| 1 | x = b + n (+ s) | `simulate_spectrum` (pipeline/chain.jl), `background`, `noise_sigma`, `line_profile` |
| 2 | smoother b̂ = S x, residual r = G_lin x on N = 24 872 bins | `SGFilter`, `smooth_unweighted`, `linear_residual` |
| 3 | iterative noise estimate: k = 0 unweighted; then S^(k+1) weighted by 1/σ̂_i^(k)², weighted windowed fit (W = 96, quadratic per window, straight line σ̂_w vs ⟨b⟩_w) or Eq. (s6); repeat until σ̂_0, σ̂_1 no longer change (tol 1e-10) | `estimate_noise`, `weighted_smooth`, `window_fit`, `lsq_noise_fit` |
| 4 | D = diag(1/σ̂_i), ρ = D r | `run_chain` |
| 5 | F r | `rfft` in `run_chain` |
| 6 | M = {f : \|G\|² ≥ τ max\|G\|²}, G(0) = 0 removed, Nyquist dropped | `retained_modes`, `retained_fraction` (κ counts ν_f) |
| 7 | divide by G(f) on M, zero elsewhere | `whitened_fourier` |
| 8 | c = D F† P_M[F r/G], z = c/√κ on the tested bins m+1…N−m (N_g = 24 744); boundary term removed with the six leading patterns of M̃ | `run_chain`, `boundary_patterns`, `remove_boundary` |
| 9 | T | `squares_T` → `RunStatistics.t_obs` |
| 10 | p from P_K = P_BC(T < t \| N_g) exp[−N_g Σ_q 2^{−(q+2)} Ψ̄_q (E_q − 1)], second order by default; Monte-Carlo p under H0^iid and H0^filt | `pk_cdf`, `pk_table`, `pvalue`, `null_iid`, `null_filt`, `mc_pvalue` |

P_BC comes from RunStatistics.jl: `squares_cdf(T, L)` for L ≤ 200 and otherwise
`squares_pvalue_approx(T, [N, L/N], ε)` with N = 100 (F(T|L) = F(T|N)^n/[1+Δ]^{n−1}). Only where
RunStatistics refuses an argument (p below its precision ε) is the same block expression
evaluated with the internal renewal recursion (`pbc_cdf_blocks`). T is always
`RunStatistics.t_obs(z, 0, 1)[1]`. If your RunStatistics version has a different call
signature, only `src/stats/squares.jl` has to be adapted.

## Installation

```
cd MadmaxRuns
julia scripts/setup.jl plots        # analysis package and plotting environment
julia --project=. -t 4 test/runtests.jl
```

## One spectrum, one SG filter

```
julia --project=. -t 32 scripts/run_pipeline.jl n=4 m=64 tau=1e-7 model=hetero level=0.05 signal=3 nnull=5000 nfilt=1000
julia --project=. -t 32 scripts/run_pipeline.jl n=4 m=64 input=myspectrum.csv out=result.jls
```

It prints every step: σ̂_0, σ̂_1 and their full iteration history, removed modes and κ, the
boundary singular values, T and its run, p from P_K, the raw residual against P_BC and the
Monte-Carlo p values. `input` is a CSV with one column x or two columns (f, x).

From Julia:

```julia
using MadmaxRuns
cfg = PipelineConfig(n = 4, m = 64, tau = 1e-7, boundary_removal = true)
S = filter_setup(cfg)                       # filter, G(f), M, κ, boundary responses, P_K table
x, b, σ, s = simulate_spectrum(Xoshiro(1), cfg; model = Heteroscedastic, level = 0.05)
r = run_chain(S, x)                         # r.noise, r.c, r.z, r.zp, r.T, r.p_pk, …
```

## Studies and figures

| Study | Paper output |
|---|---|
| `filter` | Fig. filter (fig_filter) |
| `background`, `spectra` | Figs. bkg, spectra |
| `sigmafit`, `iteration` | Figs. sigmafit (deck_nobkg), iter |
| `nobkg` | Figs. nobkg, shapes_p |
| `tauscan` | Fig. tauscan, τ_min |
| `tau0` | Figs. tau_safe_map, tau_safe_crossing, τ0 |
| `correlation` | Figs. corr, corrdiff (Δ_rms, mean over lags 1–10) |
| `fptp` | Fig. three |
| `pvalues` | Fig. pvals (four panels) |
| `examples` | Fig. examples |
| `leakage` | Fig. leakmap |
| `tolerance` | Fig. tol |
| `boundary` | Fig. bdshapes, singular values, K = 2 vs 6 check |
| `tcut` | Figs. tcut, g2nd (exact / second-order / first-order E_q against 5·10⁵ replicas) |
| `table_a`, `table_b` | Table shapes4 (formula with MC in brackets) |

```
julia --project=. -t 32 scripts/run_study.jl tauscan outdir=results
julia --project=. -t 32 scripts/run_study.jl table_b shapes=axion widths=4 outdir=results
julia --project=. -t 4  scripts/run_study.jl all quick=true outdir=smoke
julia --project=plots plots/scripts/make_figures.jl results results/figs
```

Options (`key=value`): `n m kind nspec nnull nfilt seed outdir levels models tau W filt_reestimate
taus tau_min tau0_homo tau0_hetero shapes widths amplitudes level quick`. Study-specific keywords
are passed through (e.g. `nmc=` and `Mexact=` for `tcut`, `ms=8:8:296` for `leakage`).
Defaults reproduce the paper: SG n = 4, m = 64, 5000 spectra and replicas, levels 0.12, 1, 5,
17.3 %, the 63-point τ grid, A = 0:4:28 and w = 2:2:12.

τ_min and τ0 are read by later studies from `results/tauscan.jls` and `results/tau0.jls`
(else the paper values 5.4e-5, 2.9e-4, 3.9e-4, or the values given on the command line).

## PBS

Edit `batch/config.sh` (paths, queue, `ncpus`, memory, walltime), then on the login node

```
bash batch/pbs_all.sh
```

This submits the self-test, stage 1 (12 sub-jobs), stage 2 (32 sub-jobs, the table split by
shape and width) and the final merge and figures, each depending on the previous one. Every
sub-job runs `julia -t $NCPUS` over spectra, with FFTW, BLAS and OpenMP at one thread
(`configure_threads`, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`), so the cores are not
oversubscribed. Random streams are indexed by (seed, case, spectrum), so results do not depend
on the number of threads or on how the work is split.

## Choices made in the code

- The Nyquist mode is removed together with the modes below τ (`drop_nyquist = true`), which
  gives the 97 removed modes and κ = 99.2 % of the paper at τ = 1e-7.
- H0^filt replicas go through the full chain including the iterative noise estimate
  (`filt_reestimate = true`); `false` uses the σ̂ of the data in step 4.
- Noise-only spectra (b = 0) keep σ_i of a background realisation.
- The raw residual is r/(σ̂√q) on the tested bins, tested with P_BC.
- β_max and N_g in τ0 run over the tested bins.
- In threshold scans the noise fit of a spectrum (steps 2–3, independent of τ) is done once
  and kept in full precision; steps 4–10 are repeated at every τ (`SpectrumSet`, `eval_set`).
