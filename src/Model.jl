#=
mutable struct DspModel
    # Linear constraints in CSR matrix format
    start::Vector{Cint}
    index::Vector{Cint}
    value::Vector{Float64}

    rlbd::Vector{Float64} # Row lower bounds
    rubd::Vector{Float64} # Row upper bounds

    obj::Vector{Float64}       # Linear objective cofficients
    qobjrowindex::Vector{Int}  # Quadratic objective matrix row index
    qobjcolindex::Vector{Int}  # Quadratic objective matrix column index
    qobjvalue::Vector{Float64} # Quadratic objective matrix values
    qobjnum::Int               # ??

    clbd::Vector{Float64} # Column lower bounds
    cubd::Vector{Float64} # Column upper bounds
    ctype::String         # Column types
    cname::Vector{String} # Column names

    nqrows::Int            # Number of quadratic constraints
    linnzcnt::Vector{Int}  # Number of nonzero elements in linear term for each constraint
    quadnzcnt::Vector{Int} # Number of nonzero elements in quadratic terms for each constraint
    rhs
    sense
    linstart
    linind
    linval
    quadstart
    quadrow
    quadcol
    quadval
end
=#

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

row_bounds_from_moi(set::MOI.LessThan) = (-Inf, set.upper)
row_bounds_from_moi(set::MOI.GreaterThan) = (set.lower, Inf)
row_bounds_from_moi(set::MOI.EqualTo) = (set.value, set.value)
row_bounds_from_moi(set::MOI.Interval) = @error("Interval row bounds are not supported.")
