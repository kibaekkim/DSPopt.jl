function JuMP.optimize!(m::SJ.StructuredModel; options...)

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