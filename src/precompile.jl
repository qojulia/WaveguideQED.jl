function _precompile_()
    times = 0:1:10
    bw = WaveguideBasis(1,times) 
end

using SnoopPrecompile

# precompilation causes allocation performance bugs for <v1.8 https://github.com/JuliaLang/julia/issues/35972
VERSION > v"1.8" && @precompile_setup begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        _precompile_()
    end
end