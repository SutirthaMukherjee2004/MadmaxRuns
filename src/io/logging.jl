mutable struct RunLog
    io::Union{Nothing,IOStream}
    t0::Float64
end
const LOG = RunLog(nothing, time())
open_log(path::AbstractString) = (LOG.io = open(path, "a"); LOG.t0 = time(); nothing)
close_log() = (LOG.io === nothing || close(LOG.io); LOG.io = nothing)
function log_message(msg::AbstractString)
    line = @sprintf("[%s %8.1f s] ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), time() - LOG.t0) * msg
    println(line); flush(stdout)
    LOG.io === nothing || (println(LOG.io, line); flush(LOG.io))
    return nothing
end
