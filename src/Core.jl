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