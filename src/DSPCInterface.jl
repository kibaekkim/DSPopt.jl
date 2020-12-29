"""
    Methods

List of methods available.

- `ExtensiveForm`: deterministic equivalent formulation solve
- `Benders`: Benders decomposition
- `DW`: Dantzig-Wolfe decomposition
- `Dual`: Dual decomposition
"""
@enum Methods begin
    ExtensiveForm
    Benders
    DW
    Dual
end

mutable struct DSPProblem
    p::Ptr{Cvoid}

    numRows::Dict{Int, Int}
    numCols::Dict{Int, Int}
    primVal::Float64
    dualVal::Float64
    colVal::Dict{Int, Vector{Float64}}
    rowVal::Vector{Float64}
    objective_sense::Float64
    status::Int
    solve_time::Float64

    # Number of blocks
    nblocks::Int

    # Array of block ids:
    # The size of array is not necessarily same as nblocks,
    # as block ids may be distributed to multiple processors.
    block_ids::Vector{Integer}

    is_stochastic::Bool

    dro

    # solve_type should be one of these:
    solve_type::Methods

    # MPI settings
    comm
    comm_size::Int
    comm_rank::Int

    # linear and quadratic constraints of subproblems
    is_qc::Dict{Int64,Bool}
    quadConstrs::Dict{Int64, Dict{Int64,AbstractConstraint}}
    linConstrs::Dict{Int64, Dict{Int64,AbstractConstraint}}

    function DSPProblem(fake::Bool = false)
        prob = new(
            C_NULL, # p
            Dict(), # numRows
            Dict(), # numCols
            NaN, # primVal
            NaN, # dualVal
            Dict(), # colVal
            [], # rowVal
            1., # objective_sense
            3998, # status
            0., # solve_time
            -1, # nblocks
            [], # block_ids
            false, # is_stochastic
            nothing, # dro
            DW, # solve_type
            nothing, # comm
            1, # comm_size
            0, # comm_rank
            Dict(), # is_qc
            Dict(), # quadratic constraints
            Dict() # linear constraints
        )
        if !fake
            prob.p = createEnv()
            finalizer(freeEnv, prob)
        end

        return prob
    end
end

###############################################################################
# Help functions
###############################################################################

macro dsp_ccall(func, args...)
    @static if !Sys.iswindows()
        return esc(quote
            ccall(($func, "libDsp"), $(args...))
        end)
    end
    @static if Sys.iswindows()
        return esc(quote
            ccall(($func, "libDsp"), stdcall, $(args...))
        end)
    end
end

createEnv() = @dsp_ccall("createEnv", Ptr{Cvoid}, ())

function freeEnv(dsp::DSPProblem)
    if dsp.p != C_NULL
        freeModel(dsp)
        @dsp_ccall("freeEnv", Cvoid, (Ref{Ptr{Cvoid}},), dsp.p)
        dsp.p = C_NULL
    end
    dsp.comm = nothing
    dsp.comm_size = 1
    dsp.comm_rank = 0
end

function freeModel(dsp::DSPProblem)
    @dsp_ccall("freeModel", Cvoid, (Ptr{Cvoid},), dsp.p)
    dsp.numRows = Dict()
    dsp.numCols = Dict()
    dsp.primVal = NaN
    dsp.dualVal = NaN
    dsp.colVal = Dict()
    dsp.rowVal = []
    dsp.objective_sense = 1.
    dsp.status = 3998
    dsp.solve_time = 0.
    dsp.nblocks = -1
    dsp.block_ids = []
    dsp.is_stochastic = false
    dsp.solve_type = Dual
    dsp.is_qc = Dict()
    dsp.quadConstrs = Dict()
    dsp.linConstrs = Dict()
end

###############################################################################
# Load problems
###############################################################################

readParamFile(dsp::DSPProblem, param_file::AbstractString) = @dsp_ccall(
    "readParamFile", Cvoid, 
    (Ptr{Cvoid}, Ptr{UInt8}), 
    dsp.p, param_file)

readSmps(dsp::DSPProblem, filename::AbstractString) = @dsp_ccall(
    "readSmps", Cvoid, 
    (Ptr{Cvoid}, Ptr{UInt8}), 
    dsp.p, filename)

loadFirstStage(dsp::DSPProblem, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd) = @dsp_ccall(
    "loadFirstStage", Cvoid,
    (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    dsp.p, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

loadQCQPFirstStage(dsp::DSPProblem, start, index, value, clbd, cubd, ctype, obj, qobjrowindex, qobjcolindex, qobjvalue, qobjnum, rlbd, rubd, nqrows::Int, linnzcnt::Vector{Int}, quadnzcnt::Vector{Int}, rhs::Vector{Float64}, sense::Vector{Int}, linstart::Vector{Int}, linind::Vector{Int}, linval::Vector{Float64}, quadstart::Vector{Int}, quadrow::Vector{Int}, quadcol::Vector{Int}, quadval::Vector{Float64}) = @dsp_ccall(
    "loadQCQPFirstStage", Cvoid,
    (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint},Ptr{Cint},Ptr{Cdouble}),
    dsp.p, start, index, value, clbd, cubd, ctype, obj, qobjrowindex, qobjcolindex, qobjvalue, qobjnum, rlbd, rubd, nqrows, convert(Vector{Cint}, linnzcnt), convert(Vector{Cint}, quadnzcnt), rhs, convert(Vector{Cint}, sense), convert(Vector{Cint}, linstart), convert(Vector{Cint}, linind), linval, convert(Vector{Cint}, quadstart), convert(Vector{Cint}, quadrow), convert(Vector{Cint}, quadcol), quadval)

loadSecondStage(dsp::DSPProblem, id, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd) = @dsp_ccall(
    "loadSecondStage", Cvoid,
    (Ptr{Cvoid}, Cint, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    dsp.p, id, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

loadQCQPSecondStage(dsp::DSPProblem, id, probability, start, index, value, clbd, cubd, ctype, obj, qobjrowindex, qobjcolindex, qobjvalue, qobjnum, rlbd, rubd, nqrows::Int, linnzcnt::Vector{Int}, quadnzcnt::Vector{Int}, rhs::Vector{Float64}, sense::Vector{Int}, linstart::Vector{Int}, linind::Vector{Int}, linval::Vector{Float64}, quadstart::Vector{Int}, quadrow::Vector{Int}, quadcol::Vector{Int}, quadval::Vector{Float64}) = @dsp_ccall(
    "loadQCQPSecondStage", Cvoid,
    (Ptr{Cvoid}, Cint, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint},Ptr{Cint},Ptr{Cdouble}),
    dsp.p, id, probability, start, index, value, clbd, cubd, ctype, obj, qobjrowindex, qobjcolindex, qobjvalue, qobjnum, rlbd, rubd, nqrows, convert(Vector{Cint}, linnzcnt), convert(Vector{Cint}, quadnzcnt), rhs, convert(Vector{Cint}, sense), convert(Vector{Cint}, linstart), convert(Vector{Cint}, linind), linval, convert(Vector{Cint}, quadstart), convert(Vector{Cint}, quadrow), convert(Vector{Cint}, quadcol), quadval)

loadBlockProblem(dsp::DSPProblem, id, ncols, nrows, numels, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd) = @dsp_ccall(
    "loadBlockProblem", Cvoid, (
    Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    dsp.p, id, ncols, nrows, numels, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

updateBlocks(dsp::DSPProblem) = @dsp_ccall("updateBlocks", Cvoid, (Ptr{Cvoid},), dsp.p)

###############################################################################
# solve functions
###############################################################################

for func in [:freeSolver, 
             :solveDe, 
             :solveBd, 
             :solveDd, 
             :solveDw]
    strfunc = string(func)
    @eval begin
        function $func(dsp::DSPProblem)
            return @dsp_ccall($strfunc, Cvoid, (Ptr{Cvoid},), dsp.p)
        end
    end
end

for func in [:solveBdMpi, :solveDdMpi, :solveDwMpi]
    strfunc = string(func)
    if @isdefined(MPI)
        @eval begin
            function $func(dsp::DSPProblem)
                return @dsp_ccall($strfunc, Cvoid, (Ptr{Cvoid}, MPI.MPI_Comm), dsp.p, dsp.comm)
	        end
    	end
    else
        @eval begin
            function $func(dsp::DSPProblem)
                error("MPI package is required to use this function.")
			end
		end
	end
end

###############################################################################
# Get functions
###############################################################################

for (func,rtn) in [(:getNumScenarios, Cint), 
                   (:getNumSubproblems, Cint), 
                   (:getTotalNumRows, Cint), 
                   (:getTotalNumCols, Cint), 
                   (:getStatus, Cint), 
                   (:getNumIterations, Cint), 
                   (:getNumNodes, Cint), 
                   (:getWallTime, Cdouble), 
                   (:getPrimalBound, Cdouble), 
                   (:getDualBound, Cdouble),
                   (:getNumCouplingRows, Cint),
                   (:getVersionMajor, Cint),
                   (:getVersionMinor, Cint),
                   (:getVersionPatch, Cint)]
    strfunc = string(func)
    @eval begin
        function $func(dsp::DSPProblem)
            @assert(dsp.p != C_NULL)
            return @dsp_ccall($strfunc, $rtn, (Ptr{Cvoid},), dsp.p)
        end
    end
end

function getSolution(dsp::DSPProblem, num::Integer)
    sol = zeros(num)
    if dsp.comm_rank == 0
        @dsp_ccall("getPrimalSolution", Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cdouble}), dsp.p, num, sol)
    end
    return sol
end
getSolution(dsp::DSPProblem) = getSolution(dsp, getTotalNumCols(dsp))

function getDualSolution(dsp::DSPProblem, num::Integer)
	sol = zeros(num)
    if dsp.comm_rank == 0
        @dsp_ccall("getDualSolution", Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cdouble}), dsp.p, num, sol)
    end
    return sol
end
getDualSolution(dsp::DSPProblem) = getDualSolution(dsp, getNumCouplingRows(dsp))

function getVersion(dsp::DSPProblem)
    major = getVersionMajor(dsp)
    minor = getVersionMinor(dsp)
    patch = getVersionPatch(dsp)
    return "$(major).$(minor).$(patch)"
end

###############################################################################
# Set functions
###############################################################################

setNumberOfScenarios(dsp::DSPProblem, nscen::Int) = @dsp_ccall(
    "setNumberOfScenarios", Cvoid,
    (Ptr{Cvoid}, Cint), 
    dsp.p, convert(Cint, nscen))

setDimensions(dsp::DSPProblem, ncols1::Int, nrows1::Int, ncols2::Int, nrows2::Int) = @dsp_ccall(
    "setDimensions", Cvoid,
    (Ptr{Cvoid}, Cint, Cint, Cint, Cint),
    dsp.p, convert(Cint, ncols1), convert(Cint, nrows1), convert(Cint, ncols2), convert(Cint, nrows2))
    
setIntPtrParam(dsp::DSPProblem, name::String, n::Int, v::Vector{Int}) = @dsp_ccall(
    "setIntPtrParam", Cvoid, 
    (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}),
    dsp.p, name, convert(Cint, n), convert(Vector{Cint}, v))

###############################################################################
# Other functions
###############################################################################
writeMps!(dsp::DSPProblem, filename::AbstractString) = @dsp_ccall(
    "writeMps", Cvoid, 
    (Ptr{Cvoid}, Ptr{UInt8}), 
    dsp.p, "$filename")

# end # end of module
