# Calibration of the single-spectrum mode against an ensemble of mock spectra (Sec. 2.6, 2.7, 3.7).
# For one measured spectrum, mock spectra are generated with the noise model fitted to that spectrum,
# on many background realisations, and passed through the same chain with the same settings. The
# p values of the noise-only mocks must be flat (f(p <= 0.05) ~ 0.05); the injected axion line
# shows the detection rate; the rank of the measured T among the noise-only mocks gives a p value
# that does not rely on the P_K formula. With check=true the measured residual is also compared with
# the same noise-only mocks through a set of diagnostics (normality, whiteness, noise model), which
# tests whether the measured spectrum itself obeys the assumptions the mocks are built on. Only the
# study mode generates its own mock grids; this file is used by scripts/run_pipeline.jl alone.

const CALIB_BLOCK = 128

function ask_background_count(default::Integer)
    stdin isa Base.TTY || return default
    print("How many background realisations should the spectrum be compared with? [$default]: ")
    flush(stdout)
    line = strip(readline(stdin))
    isempty(line) && return default
    n = tryparse(Int, line)
    n !== nothing && n > 0 || throw(ArgumentError("number of background realisations must be a positive integer"))
    return n
end

# Kolmogorov-Smirnov distance from the uniform distribution and its asymptotic p value.
function ks_uniform(p::AbstractVector{<:Real})
    u = sort(collect(Float64, p)); n = length(u)
    n == 0 && return (D = NaN, p = NaN)
    D = maximum(max(i / n - u[i], u[i] - (i - 1) / n) for i in 1:n)
    lambda = (sqrt(n) + 0.12 + 0.11 / sqrt(n)) * D
    q = 0.0
    for k in 1:100
        term = 2 * (-1)^(k - 1) * exp(-2 * k^2 * lambda^2)
        q += term
        abs(term) < 1e-12 && break
    end
    return (D = D, p = clamp(q, 0.0, 1.0))
end

function calibration_background_scale(grid, bhat)
    bnom = nominal_background(; grid)
    mb = median(bnom)
    isfinite(mb) && mb > 0 ||
        throw(ArgumentError("the paper background model has no positive median on this frequency grid; use bkg_model=data"))
    return median(bhat) / mb
end

function mock_noise_sigma(sigma0, sigma1, b)
    s = sigma0 .+ sigma1 .* b
    floor = 0.1 * median(abs.(s))
    floor > 0 || throw(ArgumentError("fitted noise model gives no positive noise level"))
    nfloor = count(<(floor), s)
    return max.(s, floor), nfloor
end

function calibration_line_centre(grid, m)
    off = MR.axion_offset(grid)
    b = MR.frequency_bin(grid, off)
    if !(3m < b < grid.nbins - 3m)
        b = div(grid.nbins, 2)
        off = (b - 1) * grid.df_hz
    end
    return b, off
end

# Fraction of mock chains that may fail (noise fit not converged or non-positive) before the
# calibration is declared unreliable.
const CALIB_MAX_FAIL = 0.01
# Two-sided p value below which a residual diagnostic is flagged (check mode). Each diagnostic is an
# average over ~N_g bins or ~N_g/W windows, so across the mocks it is close to Gaussian (central limit
# theorem); the flag uses the standardised deviation d = (x - mean)/std of the mocks, p = 2 Phi(-|d|).
# The empirical rank p is printed as a cross-check; it cannot go below 2/(nbkg + 1).
const CHECK_ALPHA = 1e-3
const CHECK_MIN_MOCKS = 100

# Diagnostics of one processed spectrum (check mode). z is the whitened residual after the boundary
# projection on the tested bins, rho = r / (sigma_hat sqrt(q)) the normalised raw residual. Each
# diagnostic targets one assumption of the paper; under the assumptions the noise-only mocks reproduce
# their distribution, and the measured spectrum is compared with it. All diagnostics except the two
# extreme-value notes are computed outside +-2m bins around the strongest run, in the data and in every
# mock alike, so that a narrow line - the signal searched for - does not count as a violation.
const CHECK_NAMES = ("variance", "skewness", "excess_kurtosis", "lag1_autocorr", "autocorr_rms",
                     "noise_slope", "noise_curvature", "window_dispersion", "leakage", "max_z", "neg_min_z")
const CHECK_VERDICT = ("variance", "skewness", "excess_kurtosis", "lag1_autocorr", "autocorr_rms",
                       "noise_slope", "noise_curvature", "window_dispersion", "leakage")
const CHECK_TEXT = Dict(
    "variance" => "variance of the whitened residual: the noise level is fitted correctly (Eqs. 15-18)",
    "skewness" => "skewness of the whitened residual: symmetric, Gaussian noise (Eq. 1)",
    "excess_kurtosis" => "excess kurtosis of the whitened residual: Gaussian tails, no spikes (Eq. 1)",
    "lag1_autocorr" => "lag-1 autocorrelation of the whitened residual: independent bins, filter correlation removed (Sec. 2.3)",
    "autocorr_rms" => "rms autocorrelation over lags 2-250 (the measure of Fig. 10): no longer-range correlation between bins",
    "noise_slope" => "slope of the local variance of r/sigma_hat against b_hat: noise linear in the background (Eq. 11)",
    "noise_curvature" => "curvature of the local variance of r/sigma_hat against b_hat: no non-linear noise dependence (Eq. 11)",
    "window_dispersion" => "spread of the local variance of r/sigma_hat between windows: stationary noise, no bursts",
    "leakage" => "rms of the sliding mean of r/sigma_hat over 2m+1 bins: background removed by the filter (App. B)",
    "max_z" => "largest positive excursion (a strong line also raises it)",
    "neg_min_z" => "largest negative excursion (a dip or a glitch)")
const CHECK_EPS_MAX = 0.01     # App. B: largest relative error of E[r_i^2] from the noise change in the window
const CHECK_LAGS = 250

function sliding_mean(v::AbstractVector, h::Integer)
    n = length(v); cs = [0.0; cumsum(v)]
    return [(cs[min(n, i + h) + 1] - cs[max(1, i - h)]) / (min(n, i + h) - max(1, i - h) + 1) for i in 1:n]
end

function autocorr_rms(d::AbstractVector{Float64}, lags::Integer)
    n = length(d); L = 2n
    X = MR.rfft_cached([d; zeros(L - n)])
    a = MR.irfft_cached(complex.(abs2.(X)), L)
    return sqrt(sum(abs2, view(a, 3:lags+1)) / (lags - 1)) / a[1]
end

function residual_diagnostics(S, cfg, res)
    z = res.zp; n = length(z)
    sigma = res.noise.sigma; bhat = res.noise.bhat; rows = S.rows
    rho = res.r[rows] ./ (sigma[rows] .* sqrt(S.q))
    # Exclude the strongest run and 2m bins on each side (the filter spreads a line over +-m bins), so
    # that a narrow line does not enter the noise diagnostics. The same is done in every mock.
    i0, i1 = res.run
    keep = trues(n)
    i0 > 0 && (keep[max(1, i0 - 2cfg.m):min(n, i1 + 2cfg.m)] .= false)
    nk = count(keep)
    zk = z[keep]
    mu = sum(zk) / nk
    d = (z .- mu) .* keep                       # masked bins set to zero for the correlations
    dk = zk .- mu
    m2 = sum(abs2, dk) / nk; m3 = sum(v -> v^3, dk) / nk; m4 = sum(v -> v^4, dk) / nk
    lag1 = sum(d[i] * d[i+1] for i in 1:n-1) / (sum(keep[i] & keep[i+1] for i in 1:n-1) * m2)
    W = cfg.W; nw = div(n, W)
    win = [k for k in 1:nw if all(view(keep, (k-1)*W+1:k*W))]
    s2 = [sum(abs2, view(rho, (k-1)*W+1:k*W)) / W for k in win]
    bw = [sum(view(bhat, rows[(k-1)*W+1:k*W])) / W for k in win]
    nwk = length(win)
    span = maximum(bw) - minimum(bw)
    x = span > 0 ? (bw .- sum(bw) / nwk) ./ span : zeros(nwk)
    A = [nwk sum(x) sum(abs2, x); sum(x) sum(abs2, x) sum(v -> v^3, x); sum(abs2, x) sum(v -> v^3, x) sum(v -> v^4, x)]
    coef = span > 0 ? MR.spd_solve!(Matrix{Float64}(A), [sum(s2), sum(x .* s2), sum(abs2.(x) .* s2)]) : [0.0, 0.0, 0.0]
    sc = s2 .- sum(s2) / nwk
    sm = sliding_mean(rho, cfg.m)
    leak = sqrt(sum(abs2, sm[keep]) / nk)
    return (; variance = m2, skewness = m3 / m2^1.5, excess_kurtosis = m4 / m2^2 - 3, lag1_autocorr = lag1,
            autocorr_rms = autocorr_rms(d, CHECK_LAGS) * n / nk, noise_slope = coef[2], noise_curvature = coef[3],
            window_dispersion = sqrt(sum(abs2, sc) / nwk), leakage = leak, max_z = maximum(z), neg_min_z = -minimum(z))
end

# App. B, Eq. (B.2): relative error of E[r_i^2] = q sigma_i^2 caused by the change of the fitted noise
# inside the filter window, eps_i = (1/q) sum_{0<|d|<=m} g_d^2 (sigma_(i+d)^2 / sigma_i^2 - 1), and
# Eq. (B.4)-(B.5): the largest relative change Delta_i of the fitted background across the window.
function smoothness_diagnostics(S, cfg, nf)
    m = cfg.m; g = S.F.g; q = S.q
    sfull = MR.pad_edges(nf.sigma, m); bfull = MR.pad_edges(nf.bhat, m)
    eps = 0.0; delta = 0.0
    for i in S.rows
        c = i + m
        si2 = sfull[c]^2; e = 0.0
        for k in -m:m
            k == 0 && continue
            e += g[k + m + 1]^2 * (sfull[c + k]^2 / si2 - 1)
        end
        eps = max(eps, abs(e / q))
        bc = abs(bfull[c])
        bc > 0 && (delta = max(delta, maximum(abs(bfull[c + k] - bfull[c]) for k in -m:m) / bc))
    end
    return (; eps_max = eps, delta_max = delta, passed = eps <= CHECK_EPS_MAX)
end

function empirical_two_sided(x::Real, ref::AbstractVector{<:Real})
    n = length(ref)
    hi = (1 + count(>=(x), ref)) / (n + 1); lo = (1 + count(<=(x), ref)) / (n + 1)
    return min(1.0, 2 * min(hi, lo))
end

function residual_check(data_diag, mock_diags)
    length(mock_diags) >= CHECK_MIN_MOCKS ||
        error("check needs at least $CHECK_MIN_MOCKS analysed noise-only mocks, got $(length(mock_diags)); increase nbkg")
    rows = NamedTuple[]
    for name in CHECK_NAMES
        key = Symbol(name)
        ref = [getfield(d, key) for d in mock_diags]
        x = getfield(data_diag, key)
        mu = sum(ref) / length(ref); sd = sqrt(sum(abs2, ref .- mu) / (length(ref) - 1))
        dev = sd > 0 ? (x - mu) / sd : (x == mu ? 0.0 : Inf)
        p = 2 * MR.Phi(-abs(dev))
        s = sort(ref)
        q(a) = s[clamp(round(Int, a * length(s)), 1, length(s))]
        push!(rows, (; name, value = x, mock_median = q(0.5), mock_lo = q(0.0005), mock_hi = q(0.9995), deviation = dev, p,
                       p_rank = empirical_two_sided(x, ref), verdict = name in CHECK_VERDICT, flagged = p < CHECK_ALPHA))
    end
    passed = !any(r -> r.verdict && r.flagged, rows)
    return (; rows, passed, alpha = CHECK_ALPHA)
end

# A mock chain may fail only in the ways a real noise fit can fail: no convergence, a non-positive
# noise level, or a singular weighted fit. Any other error is a bug or a resource problem and stops the
# calibration instead of being counted.
function is_fit_failure(err)
    err isa PosDefException && return true
    err isa SingularException && return true
    err isa ArgumentError || return false
    msg = err.msg
    return occursin("noise fit did not converge", msg) || occursin("fitted noise must be finite and strictly positive", msg) ||
           occursin("fitted background contains a non-finite value", msg) || occursin("noise fit has non-finite parameters", msg)
end

function calibrate_spectrum(run, S, cfg, nf, data_result; nbkg::Integer, bkg_model::Symbol, peak::Real, seed::Integer,
                            check::Bool = false, scan = nothing)
    nscan = scan === nothing ? 0 : length(scan.setups)
    T_obs = data_result.T
    grid = cfg.grid; m = cfg.m
    bmodel_full = MR.pad_edges(nf.bhat, m)
    scale = bkg_model == :paper ? calibration_background_scale(grid, nf.bhat) : 1.0
    if bkg_model == :paper && !(isapprox(grid.fref_hz, MR.PAPER_GRID.fref_hz; rtol = 1e-6) &&
                                isapprox(grid.df_hz, MR.PAPER_GRID.df_hz; rtol = 1e-6) &&
                                grid.nbins == MR.PAPER_GRID.nbins)
        log_message("calibration: the frequency grid differs from the paper grid; Eq. (27) is evaluated on this grid")
    end
    centre, offset = calibration_line_centre(grid, m)
    area = peak > 0 ? line_area_for_peak(MR.AxionLine(), 3.4, peak) : 0.0
    profile = peak > 0 ? line_profile(MR.AxionLine(), area, 3.4, grid.nbins, centre - 0.5) : zeros(grid.nbins)
    log_message(@sprintf("calibration: %d background realisations (%s), noise σ = %.6g + %.6g b, axion peak %.1fσ at f - f_ref = %.4f MHz%s",
                         nbkg, bkg_model == :paper ? "Eq. (27), r_k ~ N(0,1), scaled to the fitted median" : "fitted background of the spectrum",
                         nf.sigma0, nf.sigma1, peak, offset / 1e6, check ? ", residual check on" : ""))

    T0 = Float64[]; T1 = Float64[]; diags = NamedTuple[]; nfloored = 0
    TS0 = zeros(0, nscan); TS1 = zeros(0, nscan)
    t_start = time()
    for i0 in 1:CALIB_BLOCK:nbkg
        nblk = min(CALIB_BLOCK, nbkg - i0 + 1)
        block = MR.run_stage(run, Symbol("calibration_", i0)) do
            t0 = fill(NaN, nblk); t1 = fill(NaN, nblk); fl = zeros(Int, nblk)
            ts0 = fill(NaN, nblk, nscan); ts1 = fill(NaN, nblk, nscan)
            dg = Vector{Any}(nothing, nblk)
            Threads.@threads for j in 1:nblk
                i = i0 + j - 1
                rng = MR.stream_rng(seed + 101, i)
                b = bkg_model == :paper ? scale .* background(rng; grid) : bmodel_full
                sigma, fl[j] = mock_noise_sigma(nf.sigma0, nf.sigma1, b)
                x0 = b .+ sigma .* randn(rng, grid.nbins)
                try
                    r0 = run_chain(S, x0; compute_pvalues = false)
                    check && (dg[j] = residual_diagnostics(S, cfg, r0))
                    for k in 1:nscan
                        ts0[j, k] = run_chain(scan.setups[k], x0; noise = r0.noise, compute_pvalues = false).T
                    end
                    t0[j] = r0.T
                catch err
                    is_fit_failure(err) || rethrow()
                end
                peak > 0 || continue
                try
                    x1 = x0 .+ sigma .* profile
                    r1 = run_chain(S, x1; compute_pvalues = false)
                    for k in 1:nscan
                        ts1[j, k] = run_chain(scan.setups[k], x1; noise = r1.noise, compute_pvalues = false).T
                    end
                    t1[j] = r1.T
                catch err
                    is_fit_failure(err) || rethrow()
                end
            end
            (; T0 = t0, T1 = t1, floored = fl, diag = dg, TS0 = ts0, TS1 = ts1)
        end
        append!(T0, block.T0); append!(T1, block.T1)
        TS0 = vcat(TS0, block.TS0); TS1 = vcat(TS1, block.TS1)
        check && append!(diags, [d for d in block.diag if d !== nothing])
        nfloored += count(>(0), block.floored)
        done = i0 + nblk - 1
        (done == nbkg || div(done, CALIB_BLOCK) % 8 == 0) &&
            log_message(@sprintf("calibration: %d / %d background realisations done (%.0f s)", done, nbkg, time() - t_start))
    end
    nfail0 = count(isnan, T0); nfail1 = peak > 0 ? count(isnan, T1) : 0
    frac0 = nfail0 / nbkg; frac1 = peak > 0 ? nfail1 / nbkg : 0.0
    (nfail0 + nfail1) > 0 && log_message(@sprintf("calibration: %d noise-only and %d signal mock chains failed (noise fit not converged or not positive): %.2f %% and %.2f %% (limit %.0f %%)",
                                                 nfail0, nfail1, 100frac0, 100frac1, 100CALIB_MAX_FAIL))
    nfloored > 0 && log_message("calibration: in $nfloored mocks the fitted noise model was negative in some bins and was floored")

    ok0 = isfinite.(T0); ok1 = isfinite.(T1)
    TS0 = TS0[ok0, :]; TS1 = TS1[ok1, :]
    T0 = T0[ok0]; T1 = T1[ok1]
    isempty(T0) && error("calibration: no noise-only mock spectrum could be analysed")

    p0 = [MR.pk_pvalue(S.pk, t) for t in T0]
    cutoff = MR.pk_threshold(t -> MR.pk_cdf(S.pk, t), 0.05; lo = 1.0, hi = 150.0, tol = 1e-6)
    n0 = length(T0)
    # Failed noise-only chains are counted against the method: as false positives in the rates and as
    # realisations with T above the measured one in p_bkg. Leaving them out could thin the tail of the
    # null distribution if failures occur preferentially for extreme spectra.
    fp05 = (count(>(cutoff), T0) + nfail0) / nbkg
    fp05_ok = count(>(cutoff), T0) / n0
    fp01 = (count(<=(0.01), p0) + nfail0) / nbkg
    err05 = sqrt(0.05 * 0.95 / nbkg)
    ks = ks_uniform(p0)
    hist = MR.p_histogram(p0)
    tp05 = isempty(T1) ? NaN : count(>(cutoff), T1) / length(T1)
    n_ge = count(>=(T_obs), T0)
    p_bkg = (1 + n_ge + nfail0) / (nbkg + 1)
    p_bkg_limit = n_ge + nfail0 == 0                 # no realisation reached T: only an upper limit
    flat = abs(fp05 - 0.05) <= 3 * err05 && abs(fp05_ok - 0.05) <= 3 * err05 && ks.p > 0.01
    enough = frac0 <= CALIB_MAX_FAIL
    calibrated = flat && enough

    log_message(@sprintf("calibration: false positives f(p <= 0.05) = %.4f ± %.4f (expected 0.05; failed chains counted as positives), f(p <= 0.01) = %.4f",
                         fp05, err05, fp01))
    log_message(@sprintf("calibration: KS distance from flat p = %.4f, KS p = %.3g; p histogram (20 bins) %s",
                         ks.D, ks.p, join(hist, " ")))
    peak > 0 && log_message(@sprintf("calibration: injected axion line (peak %.1fσ) found in %.4f of the mocks%s", peak, tp05,
                                     frac1 > CALIB_MAX_FAIL ? " (more than 1 % of these chains failed; treat the rate with care)" : ""))
    log_message(p_bkg_limit ?
        @sprintf("calibration: measured T = %.4f exceeds all %d noise-only mocks: p_bkg < %.3g (Monte Carlo resolution limit, 1/(nbkg + 1))",
                 T_obs, nbkg, p_bkg) :
        @sprintf("calibration: measured T = %.4f: %d of %d noise-only mocks reach it (failed chains counted), p_bkg = %.5g",
                 T_obs, n_ge + nfail0, nbkg, p_bkg))
    if calibrated
        log_message("calibration: PASSED, the p values of this chain are flat for spectra like this one")
    else
        reasons = String[]
        flat || push!(reasons, "the noise-only p values are not flat")
        enough || push!(reasons, @sprintf("%.2f %% of the noise-only mock chains failed (limit %.0f %%)", 100frac0, 100CALIB_MAX_FAIL))
        log_message("calibration: FAILED (" * join(reasons, "; ") * "); do not trust p from P_K, use p_bkg and inspect the histogram")
    end

    scanres = nothing
    if nscan > 0
        logp(St, t) = log10(MR.tabulated_pvalue(x -> MR.pk_pvalue(St.pk, x), St.pktab, t))
        LP0 = [logp(scan.setups[k], TS0[i, k]) for i in axes(TS0, 1), k in 1:nscan]
        LP1 = [logp(scan.setups[k], TS1[i, k]) for i in axes(TS1, 1), k in 1:nscan]
        fp_tau = [(count(<=(log10(0.05)), view(LP0, :, k)) + nfail0) / nbkg for k in 1:nscan]
        tp_tau = isempty(T1) ? fill(NaN, nscan) : [count(<=(log10(0.05)), view(LP1, :, k)) / size(LP1, 1) for k in 1:nscan]
        minlp = vec(minimum(LP0; dims = 2))
        n_le = count(<=(scan.min_log10p), minlp)
        p_scan = (1 + n_le + nfail0) / (nbkg + 1)
        p_scan_limit = n_le + nfail0 == 0
        log_message("calibration: tau scan on the noise-only and signal realisations (same scan as the measured spectrum)")
        for k in 1:nscan
            log_message(@sprintf("calibration:   tau = %.1e  false positives f(p <= 0.05) = %.4f ± %.4f   injected line found in %.4f",
                                 scan.taus[k], fp_tau[k], err05, tp_tau[k]))
        end
        log_message(p_scan_limit ?
            @sprintf("calibration: tau scan: no noise-only realisation reaches the measured smallest p over the scan: p_scan < %.3g", p_scan) :
            @sprintf("calibration: tau scan: trials-corrected p_scan = %.5g (%d of %d noise-only realisations reach the measured smallest p over the scan)",
                     p_scan, n_le + nfail0, nbkg))
        scanres = (; taus = scan.taus, fp_tau, tp_tau, p_scan, p_scan_limit, min_p = 10.0^scan.min_log10p, best_tau = scan.best_tau)
    end

    checkres = nothing
    if check
        isempty(diags) && error("check: no noise-only mock could be analysed")
        dd = residual_diagnostics(S, cfg, data_result)
        checkres = residual_check(dd, diags)
        smooth = smoothness_diagnostics(S, cfg, nf)
        i0, i1 = data_result.run; width = i1 - i0 + 1
        checkres = merge(checkres, (; smooth, run_width = width, passed = checkres.passed && smooth.passed))
        log_message(@sprintf("check: measured residual against %d noise-only mocks (deviation in standard deviations of the mocks, two-sided p; flag below %.0e)",
                             length(diags), CHECK_ALPHA))
        for r in checkres.rows
            log_message(@sprintf("check:   %-18s measured %+.5g   mocks median %+.5g, 99.9%% range [%+.5g, %+.5g]   deviation %+7.2f sd, p = %.3g (rank p %.3g) %s",
                                 r.name, r.value, r.mock_median, r.mock_lo, r.mock_hi, r.deviation, r.p, r.p_rank,
                                 r.flagged ? (r.verdict ? "<- FLAGGED" : "<- note") : ""))
        end
        for r in checkres.rows
            r.flagged && log_message("check:   " * r.name * ": " * CHECK_TEXT[r.name])
        end
        log_message(@sprintf("check:   smoothness (App. B): largest relative error of the local noise variance |eps| = %.2e (limit %.0e), largest background change across the window Delta = %.3f %s",
                             smooth.eps_max, CHECK_EPS_MAX, smooth.delta_max, smooth.passed ? "" : "<- FLAGGED: the background or noise varies too fast for the filter window"))
        width > cfg.m && log_message(@sprintf("check:   note: the strongest run is %d bins wide, more than the filter half-width m = %d; the test is built for narrow lines and the filter absorbs part of a broad feature",
                                              width, cfg.m))
        log_message(checkres.passed ?
            "check: PASSED, the measured residual is consistent with the assumptions of the method (Gaussian, independent bins, noise linear in a smooth background)" :
            "check: FAILED, the measured spectrum does not follow the noise assumptions of the mocks; neither P_K nor p_bkg can be trusted for it")
    end
    return (; nbkg, bkg_model, scale, peak, line_offset_hz = offset, n_used = n0, n_failed = nfail0 + nfail1,
            n_failed_noise = nfail0, n_failed_signal = nfail1, fail_limit = CALIB_MAX_FAIL,
            n_floored = nfloored, cutoff, fp05, fp05_err = err05, fp01, ks_D = ks.D, ks_p = ks.p, hist,
            tp05, p_bkg, p_bkg_limit, calibrated, flat, T_null = sort(T0), T_signal = T1, check = checkres, scan = scanres)
end
