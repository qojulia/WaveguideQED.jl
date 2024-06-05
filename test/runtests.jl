using WaveguideQED
using Test
using SafeTestsets

function doset(descr)
    if length(ARGS) == 0
        return true
    end
    for a in ARGS
        if occursin(lowercase(a), lowercase(descr))
            return true
        end
    end
    return false
end

macro doset(descr)
    quote
        if doset($descr)
            @safetestset $descr begin include("test_"*$descr*".jl") end
        end
    end
end

println("Starting tests with $(Threads.nthreads()) threads out of `Sys.CPU_THREADS = $(Sys.CPU_THREADS)`...")

@doset "operators"
@doset "singlecavity"
# TODO doctests need fixing
#VERSION == v"1.8" && @doset "doctests"
get(ENV,"QUANTUMOPTICS_JET_TEST","")=="true" && @doset "jet"

using Aqua

doset("aqua") && begin
    Aqua.test_all(WaveguideQED,
        ambiguities=false, # TODO needs fixes
        undefined_exports=false, # TODO needs fixes
        piracy=false, # TODO needs fixes
    )
end
