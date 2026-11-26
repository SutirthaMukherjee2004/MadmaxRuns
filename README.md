# MadmaxRuns

Julia implementation of *A General Signal-Agnostic Test for Narrow Signals on Unknown Backgrounds*
(S. Mukherjee, A. Caldwell): SG background removal, noise fit, Fourier whitening, boundary-term
removal and the SQUARES statistic with the corrected distribution P_K. See the paper for the method.
The mock studies use a model of [MADMAX (MAgnetized Disc and Mirror Axion eXperiment)](https://madmax.mpp.mpg.de/news.html).

## Installation

Julia 1.12. Run all commands from the repository root.

```bash
julia scripts/setup.jl              # analysis environment
julia scripts/setup.jl plots        # + figures (CairoMakie)
julia scripts/setup.jl notebooks    # + Jupyter kernel "MadmaxRuns (Julia)"
julia --project=. -t auto test/runtests.jl
```

## Usage

### Real spectrum

```bash
julia --project=. -t auto scripts/run_pipeline.jl input=spec.csv nbkg=5000 check=true
```

Input: headerless CSV, one row per bin, `frequency_hz,power` (Hz, uniformly spaced, linear power
units). A single power column is accepted only with 25 000 bins (paper grid). Without `nbkg=` the
program asks for it in a terminal and uses 5000 otherwise. `--help` lists all options.

### Mock spectrum

```bash
julia --project=. -t auto scripts/run_pipeline.jl model=hetero level=0.05 signal=3
```

Simulates a spectrum of the paper's model (Eq. 27, Eq. 11) and runs the same chain.

### Null validation

```bash
julia --project=. -t auto scripts/validate_nulls.jl dir=path/to/signal_free_spectra
```

Runs the chain on every CSV in `dir` and reports f(p ≤ 0.05), f(p ≤ 0.01), a KS test against a flat
distribution and the p histogram; per-spectrum results go to `results/data/null_validation.csv`.

### Options

| Option | Default | |
|---|---|---|
| `input` | – | CSV spectrum; omit for a mock spectrum |
| `model`, `level`, `signal` | `homo`, `0.05`, `0` | mock spectrum: noise model, σ/b_med, axion-line peak in σ_i |
| `n`, `m`, `kind` | `4`, `64`, `SG` | SG degree, half-width, `SG` or `SGW` |
| `tau` | `1e-7` | Fourier threshold |
| `K` | `6` | boundary patterns removed (`0`: off) |
| `W`, `estimator` | `96`, `window` | noise fit window, `window` or `lsq` |
| `maxit`, `tol` | `60`, `1e-10` | noise-fit iterations and tolerance |
| `pk_order` | `2` | order of the E_q closed form |
| `seed` | `2026` | random streams |
| `calibrate` | `true` | calibration on simulated spectra |
| `nbkg` | `5000` | number of simulated spectra |
| `bkg_model` | `paper` | `paper` (Eq. 27 scaled to the fitted background) or `data` (fitted background) |
| `calib_signal` | `3` | peak of the injected axion line in σ_i; `0`: none |
| `check` | `false` | compare the residual of the spectrum with the simulated ones |
| `tauscan`, `taus` | `false`, – | repeat steps 4–10 at τ = 1e-7 … 0.2 or at a given list |
| `nnull`, `nfilt` | `0`, `0` | H0^iid and H0^filt Monte Carlo replicas |
| `output_root`, `out`, `log`, `checkpoint_dir`, `resume` | `results`, … | output locations; resume from checkpoints |

## Calibration, check and τ scan

**Calibration** (`calibrate=true`): `nbkg` spectra with backgrounds from Eq. (27) (or the fitted
background, `bkg_model=data`) and noise σ̂0 + σ̂1 b from the fit to the input, each processed with the
same settings, once without and once with an axion line.

**Check** (`check=true`): for the input and every simulated noise-only spectrum, outside ±2m bins
around the strongest run: variance, skewness, excess kurtosis, lag-1 autocorrelation and rms
autocorrelation over lags 2–250 of z; slope and curvature of the window variance of r/σ̂ against b̂;
window-to-window spread of that variance; rms of the sliding mean of r/σ̂ over 2m + 1 bins.

**τ scan** (`tauscan=true` or `taus=...`): With the
calibration on, the same scan runs on every simulated spectrum and gives the per-τ false-positive and
detection rates and p_scan, the fraction of noise-only spectra whose smallest p over the scan is at
most the measured one.

P_K is evaluated with Eq. (G.4).

## Output

Log: terminal and `results/logs/pipeline.log`.

```text
setup: SG filter, n = 4, m = 64, N = 24872, N_g = 24744, q = 0.97274, τ = 1.00e-07, 97 rfft modes removed, κ = 0.99228
steps 2–3: σ̂0 = 1.4821, σ̂1 = 0.0282591 after 4 iterations (converged: true, estimator window)
step 8: boundary singular values 57.53, 10.59, 0.5549, 0.4382, 0.1687, 0.07353
step 9: T = 23.9354 (RunStatistics.t_obs), run 3534–3539 of the tested bins = 11.540000–11.550000 MHz (6 bins)
step 10: p = 0.3049 (log10 p = -0.516) from P_K (order 2); raw residual T = 22.853, p = 0.543 from P_BC
calibration: false positives f(p <= 0.05) = 0.0520 ± 0.0069 (expected 0.05; failed chains counted as positives), f(p <= 0.01) = 0.0070
calibration: KS distance from flat p = 0.0191, KS p = 0.854; p histogram (20 bins) 52 44 52 62 39 40 58 49 58 51 57 49 50 43 49 49 45 44 57 52
calibration: injected axion line (peak 3.0σ) found in 0.9570 of the mocks
calibration: measured T = 23.9354: 292 of 1000 noise-only mocks reach it (failed chains counted), p_bkg = 0.29271
calibration: PASSED, the p values of this chain are flat for spectra like this one
check:   variance           measured +1.0043   mocks median +1.0052, 99.9% range [+0.9998, +1.0107]   deviation   -0.55 sd, p = 0.58 (rank p 0.581)
check: PASSED, ...
RESULT: CALIBRATED AND CHECKED
RESULT: T = 23.9354, strongest run 11.540000–11.550000 MHz, p(P_K) = 0.3049 (log10 -0.52)
RESULT: p_bkg = 0.2927 from 1000 background realisations
```

Tested-bin index j corresponds to row j + 2m of the input. The result is saved to
`results/data/measured_<signature>.jls`:

```julia
using MadmaxRuns
r = MadmaxRuns.load_result("results/data/measured_<signature>.jls")
r.result.T, r.result.p_pk, r.log10p, r.run_frequency_hz    # statistic, p, log10 p, candidate frequencies
r.result.noise.sigma0, r.result.noise.sigma1, r.result.zp  # noise fit, whitened residual
r.calibration                                              # fp05, ks_p, tp05, p_bkg, check, scan, ...
r.tauscan, r.summary
```

Every stage and every block of 128 simulated spectra is checkpointed under
`results/data/checkpoints/<signature>/`; rerunning the same command resumes. The signature covers
input, options and code; `resume=false` recomputes.

## Performance

32 threads:

| | before | after |
|---|---|---|
| 64 chains in parallel | 28.6 s | 1.5 s |
| calibration, 1000 spectra | ≈ 15 min | 13 s |
| real spectrum, `nbkg=5000 check=true` | ≈ 70 min | ≈ 1.5 min |
| same with `tauscan=true` | – | ≈ 6 min |

FFTW plans cached per length (`src/whitening/fftplans.jl`) instead of created per call under
FFTW's planner lock; the small normal equations of the weighted SG filter and the window fit solved
by an in-place Cholesky in Julia instead of LAPACK; T found by a linear scan of the runs (bit-identical
to `RunStatistics.t_obs`, which is still evaluated for the analysed spectrum).

## Paper studies and figures

```bash
julia --project=. scripts/run_study.jl list
julia --project=. -t 32 scripts/run_study.jl tauscan outdir=results/data resume=true
julia --project=plots plots/scripts/run_paper_notebook.jl     # figures -> results/paperplots
```

| Studies | Paper |
|---|---|
| `filter`, `background`, `spectra`, `sigmafit`, `iteration` | Figs. 1–3, 6, 7 |
| `nobkg`, `examples`, `leakage`, `tolerance`, `boundary` | Figs. 4, 5, 12, 14, B.15, F.18 |
| `tauscan`, `tau0`, `correlation`, `fptp`, `pvalues`, `tcut` | Figs. 8–11, E.16–E.17, G.19–G.20, H.21 |
| `table_a`, `table_b` | Table 1 |

The [paper figure gallery](plots/paper-reference/) contains Figures 1–21 and Table 1 as PNGs
extracted from the paper PDF.
For new results, [paperplot.ipynb](paperplot.ipynb) has one cell per figure, reads `results/data/`,
and writes generated plots to `results/paperplots/`.

The [earlier SQUARES workflow notebook](notebooks/Madmax_Runs-3.ipynb) is available separately.
Its embedded images were cleared to keep the notebook small. Run its cells to regenerate the output.

## Batch jobs

Set queue and resources in `batch/config.sh` (and `JULIA_BIN` if `julia` is not on the PATH of the
compute nodes), then:

```bash
bash batch/pbs_paper.sh       # paper studies + figures (PBS)
bash batch/pbs_all.sh         # staged studies, tables, figures (PBS)
bash batch/run_local.sh 32    # same without a scheduler
```

Logs go to `results/logs/pbs`. Stages are chained with `afterok`; `resume=true` reuses finished studies.

## Tests

```bash
julia --project=. -t auto test/runtests.jl    # 474 tests
julia --project=. -t 4 scripts/selftest.jl    # tests + quick run of every study
```

## File tree

```text
MadmaxRuns/
├── Project.toml                      analysis package: name, dependencies, compatible versions
├── Manifest.toml                     exact pinned dependency versions (Julia 1.12)
├── README.md                         this file
├── paperplot.ipynb                   notebook with one cell per paper figure (kernel "MadmaxRuns (Julia)")
├── .gitignore                        excludes results, smoke runs and checkpoints from git
├── .gitattributes                    excludes notebooks and reference figures from GitHub language statistics
├── docs/
│   └── figure-13-analysis-pipeline.png   pipeline diagram (paper Fig. 13) used by the plotting notebook
├── notebooks/
│   ├── Madmax_Runs-3.ipynb            earlier SQUARES workflow notebook
│   ├── Project.toml                  environment of the Jupyter kernel (IJulia + plots)
│   └── Manifest.toml                 its pinned versions
├── scripts/
│   ├── run_pipeline.jl               spectrum mode CLI: real (input=) or mock spectrum, full chain, RESULT
│   ├── spectrum_calibration.jl       self-calibration on simulated spectra, residual check, tau-scan calibration
│   ├── validate_nulls.jl             false-positive rate on a folder of signal-free instrument spectra
│   ├── run_study.jl                  study mode CLI: run one paper study (or all) by name
│   ├── setup.jl                      installs the environments (analysis, plots, Jupyter kernel)
│   ├── selftest.jl                   unit tests plus a quick run of every study
│   ├── merge_tables.jl               merges the per-shape/width Table 1 parts into one table
│   └── latex_table.jl                writes Table 1 as LaTeX (table_shapes4.tex)
├── src/
│   ├── MadmaxRuns.jl                 module definition: includes all files, exports the public API
│   ├── model/
│   │   ├── grid.jl                   frequency grid (25 000 bins of 2 kHz from 4.218 MHz) and bin lookup
│   │   ├── background.jl             simulated background, Eq. (27), with random parameters r_k
│   │   ├── noise.jl                  noise model sigma_i = sigma0 + sigma1 b_i, homo-/heteroscedastic, Eq. (11)
│   │   └── lineshapes.jl             Gaussian, Lorentzian, top-hat and axion (SHM) lines, Eqs. (22)-(26), (D.1)
│   ├── filter/
│   │   ├── savitzky_golay.jl         SG/SGW kernels, residual kernel g, transfer function G(f), residual r = Gx
│   │   └── weighted_sg.jl            SG filter with per-bin weights 1/sigma_i^2 (noise-fit iteration)
│   ├── noise/
│   │   └── estimate.jl               windowed and least-squares noise fits, iteration to convergence (Sec. 2.4.2)
│   ├── whitening/
│   │   ├── fftplans.jl               cached thread-safe FFT plans and the small Cholesky solver (speed)
│   │   └── fourier.jl                Fourier mask M (threshold tau), kappa, division by G(f)
│   ├── boundary/
│   │   ├── boundary.jl               boundary-term patterns of the finite spectrum and their projection (App. F)
│   │   └── tau0.jl                   lower bound tau0 on the threshold, exact and closed form (App. E)
│   ├── stats/
│   │   ├── squares.jl                SQUARES statistic T (Eq. 2) and P_BC, its distribution for i.i.d. bins (Eq. G.1)
│   │   ├── pk.jl                     P_K, the distribution corrected for removed modes (App. G), tail rules
│   │   ├── pk_exact.jl               exact E_q by Monte Carlo (Eq. G.3), for validating the closed form
│   │   ├── special.jl                chi^2 tail functions (also in log space), normal CDF and quantile
│   │   └── bivariate_normal.jl       bivariate normal CDF Phi_2 used in the second-order E_q (Eq. G.5)
│   ├── pipeline/
│   │   ├── config.jl                 PipelineConfig: all analysis settings with defaults
│   │   ├── setup.jl                  FilterSetup: everything that depends only on filter and tau (G, mask, P_K, boundary)
│   │   └── chain.jl                  run_chain (steps 2-10) and simulate_spectrum (step 1)
│   ├── nulls/
│   │   └── nulls.jl                  H0^iid and H0^filt replicas and Monte Carlo p values
│   ├── io/
│   │   ├── checkpoints.jl            stage checkpoints, run signatures, file locks (resume)
│   │   ├── results.jl                saving/loading results, CSV output, output folder layout
│   │   ├── logging.jl                time-stamped log to terminal and file
│   │   ├── cli.jl                    key=value parsing for the study CLI
│   │   └── parallel.jl               thread settings (FFTW and BLAS at one thread per Julia thread)
│   └── studies/                      one file per group of paper figures/tables (study mode)
│       ├── common.jl                 shared options, random streams, thresholds, helpers of all studies
│       ├── registry.jl               study names, quick settings, run_study with verified resume
│       ├── filter_study.jl           Fig. 1: SG kernels and transfer functions
│       ├── spectra.jl                Figs. 3, 6: backgrounds and simulated spectra
│       ├── noise_studies.jl          Figs. 2, 7: noise-fit convergence and noise estimates
│       ├── nobkg.jl                  Figs. 4, 5: SQUARES without background, four line shapes
│       ├── tables.jl                 Table 1 (a) without and (b) with background
│       ├── correlation.jl            Figs. 8-10: correlation matrices, tau scan of the autocorrelation, tau_min
│       ├── boundary_study.jl         Figs. E.16-E.17, F.18: tau0, boundary patterns and their removal
│       ├── fptp.jl                   Fig. 11: false- and true-positive rates against tau
│       ├── leakage.jl                Fig. 12: background leakage against filter degree and width
│       ├── examples.jl               Fig. 14: the full pipeline on one spectrum per noise model
│       ├── tolerance.jl              Fig. B.15: tolerance of the local-smoothness approximation (App. B)
│       ├── tcut.jl                   Figs. G.19-G.20: T distribution after the mode cut, P_K vs Monte Carlo
│       └── pvalues.jl                Fig. H.21: p-value distributions under both null hypotheses
├── plots/
│   ├── paper-reference/              figures extracted from the paper as PNGs
│   ├── Project.toml                  plotting package (CairoMakie) environment
│   ├── Manifest.toml                 its pinned versions
│   ├── src/
│   │   ├── MadmaxRunsPlots.jl        plotting module and make_all (all figures from the study results)
│   │   ├── style.jl                  paper theme, colours, labels, PDF/PNG saving
│   │   ├── figures_model.jl          filter, noise-fit convergence, background, spectra, noise estimates
│   │   ├── figures_tests.jl          no-background tests, line shapes, correlations, tau scan, FP/TP, p values
│   │   ├── figures_chain.jl          pipeline example (Fig. 14), leakage, tolerance, tau0, boundary, P_K vs Monte Carlo
│   │   └── paper.jl                  numbered paper-figure export with the required ensemble sizes
│   └── scripts/
│       ├── make_figures.jl           exports all figures as PDF and PNG
│       └── run_paper_notebook.jl     executes paperplot.ipynb and saves the numbered figures
├── batch/
│   ├── config.sh                     queue, cores, memory, walltime, paths (edit for your cluster)
│   ├── pbs_paper.sh                  submits the full paper run: studies, scans, figure export
│   ├── pbs_all.sh                    submits the staged run: self-test, stage 1, stage 2, merge and figures
│   ├── pbs_array.sh                  one PBS sub-job = one line (study + options) of a job list
│   ├── pbs_single.sh                 one PBS job running any Julia command (CMD=...)
│   ├── pbs_finish.sh                 final stage of pbs_all: table merge and figure export
│   ├── pbs_paper_finish.sh           final stage of pbs_paper: notebook and figure export
│   ├── run_local.sh                  the staged run on one machine without a scheduler
│   ├── jobs_paper_base.txt           job list: independent paper studies
│   ├── jobs_paper_scan.txt           job list: studies that need tau_min/tau0 (correlation, fptp, pvalues)
│   ├── jobs_stage1.txt               job list of stage 1 (pbs_all)
│   └── jobs_stage2.txt               job list of stage 2 (pbs_all), incl. Table 1 split by shape and width
└── test/
    ├── runtests.jl                   runs all test files below
    ├── test_model.jl                 grid, background, noise model, line shapes
    ├── test_filter.jl                SG kernels, transfer function, weighted SG vs direct fit
    ├── test_noise.jl                 noise fits and their convergence
    ├── test_whitening.jl             Fourier mask, whitening, retained fraction
    ├── test_boundary.jl              boundary patterns, projection, tau0
    ├── test_stats.jl                 T, P_BC, Phi_2, E_q moments, P_K over the full T range, fast kernels
    ├── test_chain.jl                 the full chain on simulated spectra
    ├── test_checkpoints.jl           checkpoint writing, corruption recovery, resume
    ├── test_studies.jl               quick runs of the studies, options, resume
    └── test_spectrum_mode.jl         spectrum mode end to end: calibration, check, tau scan, null validation
```

Created at run time (not in git): `results/logs/` (logs, PBS output), `results/data/` (study results,
final results, checkpoints, Table 1 CSV/LaTeX), `results/paperplots/` (figures).

## Reference

S. Mukherjee and A. Caldwell, *A General Signal-Agnostic Test for Narrow Signals on Unknown
Backgrounds*. SQUARES: F. Beaujean and A. Caldwell, J. Stat. Plan. Inference 141, 3437 (2011);
[RunStatistics.jl](https://github.com/bat/RunStatistics.jl).
