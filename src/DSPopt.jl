module DSPopt

using Pkg
using Libdl
using SparseArrays
using StructJuMP
using JuMP # To reexport, should be using (not import)
using MPI
import MathOptInterface

const SJ = StructJuMP
const MOI = MathOptInterface

@SJ.exportall JuMP
export WassersteinSet

include("DSPCInterface.jl")

function __init__()
    try
        if Libdl.find_library("libDsp") == ""
            @warn("Could not load DSP shared library. Make sure it is in your library path.")
            global dspenv = DSPProblem(true)
        else
            Libdl.dlopen("libDsp", Libdl.RTLD_GLOBAL)
            global dspenv = DSPProblem(false)
        end
    catch
        @warn("Failed to initialize the package.")
        rethrow()
    end
end

function check_dsp()
    if dspenv.p == C_NULL
        @error("DSP is not initialied.")
    end
end

include("Core.jl")
include("Model.jl")
include("Dro.jl")
include("Parallel.jl")
include("JuMP_Wrapper.jl")
include("FileIO.jl")

end # module
