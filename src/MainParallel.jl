include("./EquationDemoLibParallel.jl")
import .EquationDemoLibParallel
using Base.Threads
using Dates
function main()
    @assert(length(ARGS) == 2)
    M, N= map(x -> parse(Int, x), ARGS)
    println("Number of threads: ", nthreads())
    start_time = now()
    EquationDemoLibParallel.proximity_search(M, N)
    end_time = now()
    
    elapsed_time = Dates.value(end_time - start_time) / 1000
    println("Elapsed time: $(elapsed_time) seconds")
end

main()