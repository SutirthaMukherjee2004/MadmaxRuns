module MadmaxRuns

using LinearAlgebra, Random, Statistics, Printf, Serialization, DelimitedFiles, Dates
using FFTW
using SHA, CodecZlib
using SpecialFunctions: erf, erfc, erfcinv, loggamma, gamma_inc
import RunStatistics

include("model/grid.jl")
include("model/background.jl")
include("model/noise.jl")
include("model/lineshapes.jl")

include("whitening/fftplans.jl")
include("filter/savitzky_golay.jl")
include("filter/weighted_sg.jl")

include("noise/estimate.jl")

include("whitening/fourier.jl")

include("boundary/boundary.jl")
include("boundary/tau0.jl")

include("stats/special.jl")
include("stats/bivariate_normal.jl")
include("stats/squares.jl")
include("stats/pk.jl")
include("stats/pk_exact.jl")

include("io/logging.jl")
include("io/results.jl")
include("io/checkpoints.jl")

include("pipeline/config.jl")
include("pipeline/setup.jl")
include("pipeline/chain.jl")

include("nulls/nulls.jl")

include("io/parallel.jl")

include("studies/common.jl")
include("studies/filter_study.jl")
include("studies/nobkg.jl")
include("studies/tables.jl")
include("studies/spectra.jl")
include("studies/noise_studies.jl")
include("studies/correlation.jl")
include("studies/fptp.jl")
include("studies/pvalues.jl")
include("studies/examples.jl")
include("studies/leakage.jl")
include("studies/tolerance.jl")
include("studies/boundary_study.jl")
include("studies/tcut.jl")
include("studies/registry.jl")
include("io/cli.jl")

export SpectrumGrid, PAPER_GRID, frequencies, background, nominal_background, median_background,
       NoiseModel, Homoscedastic, Heteroscedastic, noise_parameters, noise_sigma, parse_noise_model,
       LineShape, GaussianLine, LorentzianLine, TopHatLine, AxionLine, ALL_SHAPES, line_cdf,
       line_profile, line_area_for_peak, shape_name, frequency_bin
export SGFilter, sg_kernel, sgw_kernel, residual_kernel, transfer_function, linear_residual,
       linear_residual!, smooth_unweighted, weighted_smooth, weighted_smooth_direct, kernel_moments,
       noise_pass_fraction
export NoiseFit, estimate_noise, window_fit, lsq_noise_fit, sigma_profile
export retained_modes, removed_modes, retained_fraction, whiten, whiten!, whitened_fourier
export BoundaryBasis, boundary_responses, boundary_patterns, remove_boundary, boundary_end_bins,
       tau0_exact, tau0_closed, boundary_significance
export squares_T, squares_T_run, pbc_cdf, pbc_pvalue, pbc_cdf_recursion, gamma_lags, correlation_matrix,
       Eq_first, Eq_second, Eq_exact, pk_cdf, pk_pvalue, pk_log10pvalue, PKTable, pk_table, pvalue, pk_threshold,
       Psi, Psibar, Phi, Phiinv, Phi2
export PipelineConfig, FilterSetup, filter_setup, ChainResult, run_chain, analyse_spectrum,
       simulate_spectrum, inject_line!
export NullSet, null_iid, null_filt, mc_pvalue
export run_study, STUDIES, STUDY_ORDER, StudyOptions, save_result, load_result, log_message, open_log, close_log,
       configure_threads, parse_cli, merge_tables, tau_setups, SpectrumSet, fit_set, eval_set,
       prepare_output_layout, with_run_log

end
