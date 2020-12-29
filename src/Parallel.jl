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
myblocks(nblocks::Int) = getBlockIds(nblocks)

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

"""
    parallelize

Initialize DSP environment related to MPI, if MPI is initialized.
"""
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

"""
    myrank

Returns the MPI rank (i.e., process id), starting from 0
"""
myrank() = dspenv.comm_rank

"""
    mysize

Returns the number of MPI processes
"""
mysize() = dspenv.comm_size