function configure_threads(; fftw_threads::Integer = 1, blas_threads::Integer = 1)
    FFTW.set_num_threads(fftw_threads)
    BLAS.set_num_threads(blas_threads)
    log_message("Julia threads: $(Threads.nthreads()), FFTW threads: $fftw_threads, BLAS threads: $blas_threads")
    return nothing
end

function tmap(f, n::Integer)
    out = Vector{Any}(undef, n)
    Threads.@threads for i in 1:n
        out[i] = f(i)
    end
    return [o for o in out]
end
