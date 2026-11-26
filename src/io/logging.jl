mutable struct RunLog
    io::Union{Nothing,IOStream}
    t0::Float64
end
const LOG = RunLog(nothing, time())
function open_log(path::AbstractString)
    mkpath(dirname(abspath(path)))
    close_log()
    LOG.io = open(path, "a")
    LOG.t0 = time()
    return nothing
end
close_log() = (LOG.io === nothing || close(LOG.io); LOG.io = nothing)
function log_message(msg::AbstractString)
    line = @sprintf("[%s %8.1f s] ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), time() - LOG.t0) * msg
    println(line); flush(stdout)
    LOG.io === nothing || (println(LOG.io, line); flush(LOG.io))
    return nothing
end

function with_run_log(f, path::AbstractString)
    mkpath(dirname(abspath(path)))
    previous = (; io = LOG.io, t0 = LOG.t0)
    return open(path, "a") do io
        LOG.io = io
        LOG.t0 = time()
        try
            f()
        catch err
            log_message("Run stopped: $(sprint(showerror, err))")
            rethrow()
        finally
            LOG.io = previous.io
            LOG.t0 = previous.t0
        end
    end
end
