module DSPopt

@enum Methods begin
    ExtensiveForm
    Benders
    DW
    Dual
end

using Libdl
using StructJuMP
using SparseArrays
import MathOptInterface

const SJ = StructJuMP
const MOI = MathOptInterface

using JuMP # To reexport, should be using (not import)
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

include("Dro.jl")

function SJ.optimize!(m::SJ.StructuredModel; options...)

    # free any existing model pointer
    freeModel(dspenv)

    # set options
    setoptions!(options)

    # classify the second stage constraints into linear ones and quadratic ones
    classifyConstrs(m)

    # load problem
    load_problem!(m)
    
    # solve DSP problem
    solve!()
    post_solve!()
    
    return termination_status(m)
end

" write model into an MPS file "
function writeMps!(m::SJ.StructuredModel, filename::AbstractString; options...)
    # free any existing model pointer
    freeModel(dspenv)

    # set options
    setoptions!(options)

    # classify the second stage constraints into linear ones and quadratic ones
    classifyConstrs(m)

    # load problem
    load_problem!(m)

    # write problem into MPS file
    writeMps!(dspenv, filename)

    return
end

function classifyConstrs(m::SJ.StructuredModel)
    
    nchildren = length(SJ.getchildren(m))
    dspenv.is_qc = Dict{Int64,Bool}()
    
    # first-stage quadratic constraints
    dspenv.is_qc[0] = false

    numQc = 0
    numLc = 0

    dspenv.quadConstrs[0] = Dict{Int64,AbstractConstraint}()
    dspenv.linConstrs[0] = Dict{Int64,AbstractConstraint}()
    quadConstrs_temp = Dict{Int64,AbstractConstraint}()
        
    # partition the constraint set into linear constraints and quadratic constraints
    ## count the number of linear and quadratic constraints and adjust the indices of linear constraints so that its keys = 1:numLc[id]
    for i = 1:length(m.constraints)
        c = m.constraints[i]
        if typeof(c.func) <: GenericAffExpr
            numLc += 1
            dspenv.linConstrs[0][numLc] = c
        elseif (typeof(c.func) <: GenericQuadExpr)
            numQc += 1
            quadConstrs_temp[numQc] = c
        else
            @error("Current version accepts constraints that are either quadratic or affine.")
        end
    end
        
    if numQc > 0
        dspenv.is_qc[0] = true
        # adjust the indices of quadratic constraints so that its keys = numLc+1:nrows_core
        for i = 1:numQc
            dspenv.quadConstrs[0][i + numLc] = quadConstrs_temp[i]
        end

        m.constraints = merge(dspenv.linConstrs[0], dspenv.quadConstrs[0])
        print("Core numLc: ", numLc, ", numQc: ", numQc, "\n")
    end

    # second-stage quadratic quadratic constraints
    for (id, blk) in SJ.getchildren(m)

        dspenv.is_qc[id] = false
        
        numQc = 0
        numLc = 0
        
        dspenv.quadConstrs[id] = Dict{Int64,AbstractConstraint}()
        dspenv.linConstrs[id] = Dict{Int64,AbstractConstraint}()
        quadConstrs_temp = Dict{Int64,AbstractConstraint}()
        
        # partition the constraint set into linear constraints and quadratic constraints
        ## count the number of linear and quadratic constraints and adjust the indices of linear constraints so that its keys = 1:numLc[id]
        for i = 1:length(blk.constraints)
            c = blk.constraints[i]
            if typeof(c.func) <: GenericAffExpr
                numLc += 1
                dspenv.linConstrs[id][numLc] = c
            elseif (typeof(c.func) <: GenericQuadExpr)
                numQc += 1
                quadConstrs_temp[numQc] = c
            else
                @error("Current version accepts constraints that are either quadratic or affine.")
            end
        end
        
        if numQc > 0
            dspenv.is_qc[id] = true
            # adjust the indices of quadratic constraints so that its keys = numLc+1:nrows_core
            for i = 1:numQc
                dspenv.quadConstrs[id][i + numLc] = quadConstrs_temp[i]
            end

            blk.constraints = merge(dspenv.linConstrs[id], dspenv.quadConstrs[id])
            # print("Scen ", id, " numLc: ", numLc, ", numQc: ", numQc, "\n")
        end
    end
end

function get_model_data(m::SJ.StructuredModel, id::Int = 0)

    # retrieve classified constraints
    linConstr = dspenv.linConstrs[id]
    quadConstrs = dspenv.quadConstrs[id]

    # Get a column-wise sparse matrix
    start, index, value, rlbd, rubd = get_constraint_matrix(m, linConstr)

    # column information
    clbd = Vector{Float64}(undef, num_variables(m))
    cubd = Vector{Float64}(undef, num_variables(m))
    ctype = ""
    cname = Vector{String}(undef, num_variables(m))
    for i in 1:num_variables(m)
        vref = SJ.StructuredVariableRef(m, i)
        v = m.variables[vref.idx]
        if v.info.integer
            ctype = ctype * "I"
        elseif v.info.binary
            ctype = ctype * "B"
        else
            ctype = ctype * "C"
        end
        if v.info.has_fix
            clbd[vref.idx] = v.info.fixed_value
            cubd[vref.idx] = v.info.fixed_value
        elseif v.info.binary
            clbd[vref.idx] = 0.0
            cubd[vref.idx] = 1.0
        else
            clbd[vref.idx] = v.info.has_lb ? v.info.lower_bound : -Inf
            cubd[vref.idx] = v.info.has_ub ? v.info.upper_bound : Inf
        end
        cname[vref.idx] = m.varnames[vref.idx]
    end

    # objective coefficients
    obj = zeros(num_variables(m))
    if !(objective_function_type(m) <: Real)
        for (v,coef) in objective_function(m).terms
            obj[v.idx] = coef
        end
    end

    if objective_sense(m) == MOI.MAX_SENSE
        dspenv.objective_sense = -1.
        obj .*= -1
    end

    nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval = get_qc_data(m, quadConstrs)
    return start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval
end

function get_constraint_matrix(m::SJ.StructuredModel, linConstrs::Dict{Int64,AbstractConstraint})

    is_parent = SJ.getparent(m) == nothing ? true : false

    num_rows = 0 # need to count
    num_cols = num_variables(m)
    if !is_parent
        num_cols += num_variables(SJ.getparent(m))
    end
    
    # count the number of nonzero elements
    nnz = 0
    for (i,cons) in linConstrs
        nnz += length(cons.func.terms)
        num_rows += 1
    end
    
    rind = Vector{Int}(undef, nnz)
    cind = Vector{Int}(undef, nnz)
    value = Vector{Float64}(undef, nnz)
    rlbd = Vector{Float64}(undef, num_rows)
    rubd = Vector{Float64}(undef, num_rows)
    
    pos = 1
    for (i,cons) in linConstrs
        for (var,coef) in cons.func.terms
            rind[pos] = i
            if is_parent
                cind[pos] = var.idx
            elseif JuMP.owner_model(var) == SJ.getparent(m)
                cind[pos] = var.idx
            else
                cind[pos] = var.idx + num_variables(SJ.getparent(m))
            end
            value[pos] = coef
            pos += 1
        end

        # row bounds
        rlbd[i], rubd[i] = row_bounds_from_moi(cons.set)
    end
    @assert(pos-1==nnz)
    
    # NOTE: DSP takes CSR (row-wise sparse matrix) format.
    # So I simply switch the entries for columns and rows.
    mat = sparse(cind, rind, value, num_cols, num_rows)
    dropzeros!(mat)

    # sparse description
    start = convert(Vector{Cint}, mat.colptr .- 1)
    index = convert(Vector{Cint}, mat.rowval .- 1)
    value = mat.nzval

    return start, index, value, rlbd, rubd
end

row_bounds_from_moi(set::MOI.LessThan) = (-Inf, set.upper)
row_bounds_from_moi(set::MOI.GreaterThan) = (set.lower, Inf)
row_bounds_from_moi(set::MOI.EqualTo) = (set.value, set.value)
row_bounds_from_moi(set::MOI.Interval) = @error("Interval row bounds are not supported.")

function get_qc_data(m::SJ.StructuredModel, quadConstrs::Dict{Int64,AbstractConstraint})
    
    is_parent = SJ.getparent(m) == nothing ? true : false

    nqrows = length(quadConstrs)
    
    linnzcnt = Vector{Int}(undef, nqrows)
    quadnzcnt = Vector{Int}(undef, nqrows)

    linstart = Vector{Int}(undef, nqrows)
    quadstart = Vector{Int}(undef, nqrows)

    rhs = Vector{Float64}(undef, nqrows)
    sense = Vector{Int}(undef, nqrows)

    i = 1

    for (key,con) in quadConstrs
        
        # Remove terms with `0` coefficients.
        drop_zeros!(con.func)
        
        # affine terms
        linnzcnt[i] = length(con.func.aff.terms)
        
        if i > 1
            linstart[i] = linstart[i-1] + linnzcnt[i-1]
        else 
            linstart[i] = 0
        end
        
        # print(i, "th linnzcnt: ", linnzcnt[i], "\n")

        # for (var,coef) in con.func.aff.terms
        #     print("lin var: " , var, ", coef: ", coef, "\n")
        # end

        # quad terms
        quadnzcnt[i] = length(con.func.terms);

        if i > 1
            quadstart[i] = quadstart[i-1] + quadnzcnt[i-1]
        else 
            quadstart[i] = 0
        end

        # print(i, "th quadnzcnt: ", quadnzcnt[i], "\n")
        
        # for (var,coef) in con.func.terms
        #     print("quad col: " , var.a, "\n")
        #     print("quad row: " ,var.b, "\n")
        #     print("coef: ", coef, "\n")
        # end

        if typeof(con.set) <: MOI.LessThan
            sense[i] = 'L'
            # sense = sense * "L"
            rhs[i] = con.set.upper
        elseif typeof(con.set) <: MOI.GreaterThan
            sense[i] = 'G'
            # sense = sense * "G"
            rhs[i] = con.set.lower
        elseif typeof(con.set) <: MOI.EqualTo
            @error("quadratic constraints are either 'LessThan' or 'GreaterThan'.")
        elseif typeof(con.set) <: MOI.Interval
            @error("quadratic constraints are either 'LessThan' or 'GreaterThan'.")
        end
        
        i += 1
    end

    @assert(i == nqrows + 1)

    total_linnzcnt = nqrows == 0 ? 0 : linstart[nqrows] + linnzcnt[nqrows]
    total_quadnzcnt = nqrows == 0 ? 0 : quadstart[nqrows] + quadnzcnt[nqrows]

    linind = Vector{Int}(undef, total_linnzcnt)
    linval = Vector{Float64}(undef, total_linnzcnt)

    quadrow = Vector{Int}(undef, total_quadnzcnt)
    quadcol = Vector{Int}(undef, total_quadnzcnt)
    quadval = Vector{Float64}(undef, total_quadnzcnt)

    linpos = 1
    quadpos = 1
    for (key,con) in quadConstrs
        # affine terms
        for (var,coef) in con.func.aff.terms
            if is_parent
                linind[linpos] = var.idx - 1
            else 
                if JuMP.owner_model(var) == SJ.getparent(m)
                    linind[linpos] = var.idx - 1
                else 
                    linind[linpos] = var.idx - 1 + num_variables(SJ.getparent(m))
                end
            end
            linval[linpos] = coef
            linpos += 1
            # print("lin var: " , var, ", coef: ", coef, "\n")
        end
        # quad terms
        for (var,coef) in con.func.terms
            if is_parent
                quadrow[quadpos] = var.a.idx - 1 
                quadcol[quadpos] = var.b.idx - 1 
            else 
                if JuMP.owner_model(var.a) == SJ.getparent(m)
                    quadrow[quadpos] = var.a.idx - 1 
                else 
                    quadrow[quadpos] = var.a.idx - 1 + num_variables(SJ.getparent(m))
                end
                if JuMP.owner_model(var.b) == SJ.getparent(m)
                    quadcol[quadpos] = var.b.idx - 1 
                else 
                    quadcol[quadpos] = var.b.idx - 1 + num_variables(SJ.getparent(m))
                end
            end
            quadval[quadpos] = coef
            quadpos += 1
            # print("quad col: " , var.a, "\n")
            # print("quad row: " ,var.b, "\n")
            # print("coef: ", coef, "\n")
        end
    end

    if nqrows > 0
        @assert(linpos-1==sum(linnzcnt[k] for k=1:nqrows))
        @assert(quadpos-1==sum(quadnzcnt[k] for k=1:nqrows))
    end

    return nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval
end

"""
    setoptions!

This takes `options` arguments for setting DSP options.

# Arguments
- `options`: possible keys are `:param` (path to parameter file), 
    `:is_stochastic` (`true` if stochastic model; otherwise, `false`), 
    and `:solve_type` (algorithm type; see `Methods`).
"""
function setoptions!(options)
    check_dsp()
    for (optname, optval) in options
        if optname == :param
            readParamFile(dspenv, optval)
        elseif optname == :is_stochastic
            dspenv.is_stochastic = optval
        elseif optname == :solve_type
            if optval in instances(Methods)
                dspenv.solve_type = optval
            else
                @warn("solve_type $optval is not available.")
            end
        else
            @warn("Options $optname is not available.")
        end
    end
end

"""
    load_problem!

Load problem from StructJuMP

# Arguments
- `m`: StructJuMP model
"""
function load_problem!(m::SJ.StructuredModel)

    # set number of blocks (scenarios)
    dspenv.nblocks = SJ.num_scenarios(m)

    if dspenv.is_stochastic
        loadStochasticProblem!(m)
    else
        loadStructuredProblem!(m)
    end

    setBlocks()
end

"""
    loadStochasticProblem!

Load stochastic programming problem from StructJuMP

# Arguments
- `model`: StructJuMP model
"""
function loadStochasticProblem!(model::SJ.StructuredModel)

    nscen = dspenv.nblocks
    ncols1 = length(model.variables)
    nrows1 = length(dspenv.linConstrs[0])
    ncols2 = 0
    nrows2 = 0
    for (id, subm) in SJ.getchildren(model)
        ncols2 = length(subm.variables)
        nrows2 = length(dspenv.linConstrs[id])
        break
    end

    # set scenario indices for each MPI processor
    if dspenv.comm_size > 1
        ncols2 = MPI.Allreduce([ncols2], MPI.MAX, dspenv.comm)[1]
        nrows2 = MPI.Allreduce([nrows2], MPI.MAX, dspenv.comm)[1]
    end

    # set DSPProblem data
    dspenv.numCols[0] = ncols1
    dspenv.numRows[0] = nrows1
    dspenv.colVal[0] = Vector{Float64}(undef, ncols1)
    for (id, blk) in SJ.getchildren(model)
        dspenv.numCols[id] = ncols2
        dspenv.numRows[id] = nrows2
        dspenv.colVal[id] = Vector{Float64}(undef, ncols2)
    end

    setNumberOfScenarios(dspenv, nscen)
    setDimensions(dspenv, ncols1, nrows1, ncols2, nrows2)

    qc_supported = true

    # set problem data
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval = get_model_data(model)
    if nqrows == 0
        loadFirstStage(dspenv, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
    else
        if getVersionMajor(dspenv) >= 2
            loadQCQPFirstStage(dspenv, start, index, value, clbd, cubd, ctype, obj, C_NULL, C_NULL, C_NULL, 0, rlbd, rubd, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval)
        else
            loadFirstStage(dspenv, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
            qc_supported = false
        end
    end
    
    for (id, blk) in SJ.getchildren(model)
        probability = SJ.getprobability(model)[id]
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval = get_model_data(blk, id)
        if nqrows == 0
            loadSecondStage(dspenv, id-1, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
        else
            if getVersionMajor(dspenv) >= 2
                loadQCQPSecondStage(dspenv, id-1, probability, start, index, value, clbd, cubd, ctype, obj, C_NULL, C_NULL, C_NULL, 0, rlbd, rubd, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval)
            else
                loadSecondStage(dspenv, id-1, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
                qc_supported = false
            end
        end
    end

    if qc_supported == false
        @warn "QCQP is not supported with DSP version $(getVersion(dspenv)). The unsupported objective/constraints are ignored."
    end

    # Set DRO data
    loadDroData(dspenv, dspenv.dro)
end

"""
    loadStructuredProblem!

Load generic block-structured problem from StructJuMP

# Arguments
- `model`: StructJuMP model
"""
function loadStructuredProblem!(model::SJ.StructuredModel)

    ncols1 = length(model.variables)
    nrows1 = length(model.constraints)

    # set DSPProblem data
    dspenv.numCols[0] = ncols1
    dspenv.numRows[0] = nrows1
    dspenv.colVal[0] = Vector{Float64}(undef, ncols1)
    for (id, blk) in SJ.getchildren(model)
        ncols2 = length(blk.variables)
        nrows2 = length(dspenv.linConstrs[id])
        dspenv.numCols[id] = ncols2
        dspenv.numRows[id] = nrows2
        dspenv.colVal[id] = Vector{Float64}(undef, ncols2)
    end
    
    # load master
    start, index, value, rlbd, rubd, obj1, clbd1, cubd1, ctype1, cname, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval = get_model_data(model)
    loadBlockProblem(dspenv, 0, ncols1, nrows1, start[end],
        start, index, value, clbd1, cubd1, ctype1, obj1, rlbd, rubd)

    # Check if the master block has quadratic constraints.
    has_qc = nqrows == 0 ? false : true

    # going over blocks
    for (id, blk) in SJ.getchildren(model)
        ncols2 = length(blk.variables)
        nrows2 = length(dspenv.linConstrs[id])
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval = get_model_data(blk, id)
        loadBlockProblem(dspenv, id, ncols1 + ncols2, nrows2, start[end], 
            start, index, value, [clbd1; clbd], [cubd1; cubd], [ctype1; ctype], [obj1; obj], rlbd, rubd)

        # Check if the sub-block has quadratic constraints.
        if !has_qc && nqrows > 0
            has_qc = true
        end
    end

    if has_qc
        @warn "Quadratic constraints are not supported for generic structured programs and thus will be ignored."
    end

    # Finalize loading blocks
    updateBlocks(dspenv)
end

function solve!()
    if dspenv.comm_size == 1
        if dspenv.solve_type == Dual
            if dspenv.is_stochastic
                solveDd(dspenv);
            else
                @warn("Dual decomposition is available for stochastic programming only.")
                return
            end
        elseif dspenv.solve_type == Benders
            if dspenv.is_stochastic
                solveBd(dspenv);
            else
                @warn("Benders decomposition is available for stochastic programming only.")
                return
            end
        elseif dspenv.solve_type == ExtensiveForm
            solveDe(dspenv);
        elseif dspenv.solve_type == DW
            solveDw(dspenv);
        else
            @error("Unexpected error")
            return
        end
    elseif dspenv.comm_size > 1
        if dspenv.solve_type == Dual
            if dspenv.is_stochastic
                solveDdMpi(dspenv);
            else
                @warn("Dual decomposition is available for stochastic programming only.")
                return
            end
        elseif dspenv.solve_type == Benders
            if dspenv.is_stochastic
                solveBdMpi(dspenv);
            else
                @warn("Benders decomposition is available for stochastic programming only.")
                return
            end
        elseif dspenv.solve_type == DW
            solveDwMpi(dspenv);
        elseif dspenv.solve_type == ExtensiveForm
            solveDe(dspenv);
        else
            @error("Unexpected error")
            return
        end
    end

    # solution status
    dspenv.status = getStatus(dspenv)
end

function post_solve!()
    if dspenv.status == 3998
        return
    end

    # get solution time
    dspenv.solve_time = getWallTime(dspenv)

    # primal and dual objective values
    dspenv.primVal = getPrimalBound(dspenv) * dspenv.objective_sense
    dspenv.dualVal = getDualBound(dspenv) * dspenv.objective_sense

    if dspenv.solve_type == Dual
        dspenv.rowVal = getDualSolution(dspenv)
    end

    primVal = dspenv.primVal
    if mysize() > 1
        primVal = MPI.bcast(dspenv.primVal, 0, dspenv.comm)
    end

    if abs(primVal) < 1.0e+20
        primsol = getSolution(dspenv)

        # parse solution to each block
        n_start = 1
        n_end = dspenv.numCols[0]
        dspenv.colVal[0] = primsol[n_start:n_end]
        n_start += dspenv.numCols[0]

        numBlockCols = getNumBlockCols()
        for i in 1:dspenv.nblocks
            n_end += numBlockCols[i]
            dspenv.colVal[i] = primsol[n_start:n_end]
            n_start += numBlockCols[i]
        end
    end
end

"""
    getBlockIds

Get the vector of block IDs that are assigned to the current MPI rank

# Arguments
- `nblocks`: number of blocks (by default, `DSPProblem.nblocks`)
- `master_has_subblocks`: indicate whether the master process solves (scenario) blocks or not (should always be `false`). This should not be modified.
"""
function getBlockIds(nblocks::Int = dspenv.nblocks; master_has_subblocks::Bool = true)::Vector{Int}

    if getVersionMajor(dspenv) >= 2
        master_has_subblocks = false
    end

    # processor info
    mysize = dspenv.comm_size
    myrank = dspenv.comm_rank
    # empty block ids
    proc_idx_set = Int[]
    # DSP is further parallelized with mysize > nblocks.
    modrank = myrank % nblocks
    # If we have more than one processor,
    # do not assign a sub-block to the master.
    if master_has_subblocks
        # assign sub-blocks in round-robin fashion
        for s = modrank:mysize:(nblocks-1)
            push!(proc_idx_set, s+1)
        end
    else
        if mysize > 1
            if myrank == 0
                return proc_idx_set
            end
            # exclude master
            mysize -= 1;
            modrank = (myrank-1) % nblocks
        end
        # assign sub-blocks in round-robin fashion
        for s = modrank:mysize:(nblocks-1)
            push!(proc_idx_set, s+1)
        end
    end

    # return assigned block ids
    return proc_idx_set
end

"""
    getNumBlockCols

Get the number of columns for each block
"""
function getNumBlockCols()::Dict{Int,Int}
    # get number of block columns
    numBlockCols = Dict{Int,Int}()
    if dspenv.comm_size > 1
        num_proc_blocks = convert(Vector{Cint}, MPI.Allgather(length(dspenv.block_ids), dspenv.comm))
        #@show num_proc_blocks
        #@show collect(keys(blocks))
        block_ids = MPI.Allgatherv(convert(Vector{Cint}, dspenv.block_ids), num_proc_blocks, dspenv.comm)
        #@show block_ids
        ncols_to_send = Cint[dspenv.numCols[id] for id in dspenv.block_ids]
        #@show ncols_to_send
        ncols = MPI.Allgatherv(ncols_to_send, num_proc_blocks, dspenv.comm)
        #@show ncols
        @assert(length(block_ids) == dspenv.nblocks)
        for i in 1:dspenv.nblocks
            numBlockCols[block_ids[i]] = ncols[i]
        end
    else
        for id in dspenv.block_ids
            numBlockCols[id] = dspenv.numCols[id]
        end
    end
    return numBlockCols
end

"""
    setBlocks

The function sets the number of blocks (e.g., scenarios for stochastic program) and their Ids.
"""
function setBlocks()
    dspenv.block_ids = getBlockIds()
    @dsp_ccall("setIntPtrParam", Cvoid, (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}),
        dspenv.p, "ARR_PROC_IDX", convert(Cint, length(dspenv.block_ids)), convert(Vector{Cint}, dspenv.block_ids .- 1))
end

function parallelize(comm)
    if @isdefined(MPI)
        if MPI.Initialized()
            dspenv.comm = comm
            dspenv.comm_size = MPI.Comm_size(dspenv.comm)
            dspenv.comm_rank = MPI.Comm_rank(dspenv.comm)
        else
            @error("MPI is not initialized.")
        end
    else
        @error("MPI.jl is not defined.")
    end
end

myrank() = dspenv.comm_rank
mysize() = dspenv.comm_size
myblocks(nblocks::Int) = getBlockIds(nblocks)

#########################
# Redefine JuMP functions
#########################

function JuMP.termination_status(m::StructuredModel)
    status = MOI.OPTIMIZE_NOT_CALLED
    if dspenv.status == 3000
        status = MOI.OPTIMAL
    elseif dspenv.status == 3005
        status = MOI.ALMOST_OPTIMAL
    elseif dspenv.status == 3001
        status = MOI.INFEASIBLE
    elseif dspenv.status == 3002
        status = MOI.DUAL_INFEASIBLE
    elseif dspenv.status in [3004,3007]
        status = MOI.TIME_LIMIT
    elseif dspenv.status == 3010
        status = MOI.ITERATION_LIMIT
    elseif dspenv.status in [3014,3015]
        status = MOI.OBJECTIVE_LIMIT
    elseif dspenv.status == 3006
        status = MOI.NODE_LIMIT
    elseif dspenv.status in [3008,3009,3016]
        status = MOI.OTHER_LIMIT
    else
        status = MOI.OTHER_ERROR
    end
    return status
end

JuMP.objective_value(m::StructuredModel) = dspenv.primVal
JuMP.dual_objective_value(m::StructuredModel) = dspenv.dualVal

function JuMP.value(v::SJ.StructuredVariableRef; result = 1)::Float64
    if !isnothing(v.model.parent)
        @warn("Solution values are available for the parent model only. Please use value() to access all children solutions in vector form.")
        return NaN
    else
        return dspenv.colVal[0][v.idx]
    end
end
JuMP.value() = dspenv.colVal
JuMP.dual() = dspenv.rowVal

JuMP.solve_time(m::StructuredModel) = dspenv.solve_time

end # module
