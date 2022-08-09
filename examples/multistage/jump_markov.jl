using JuMP
using GLPK

T = 3; # number of stages

# scenario lattice w/ T stages
nodes = 1:2T-1; # nodes in lattice
stage = Vector{Int64}(undef, 2T-1); # stage of each node in lattice
d = Vector{Int64}(undef, 2T-1); # demand vector
for i in nodes
    if i == 1 # stage/demand for root node
        stage[i] = 0;
        d[i] = 1;
    else # stage/demand for all other nodes
        stage[i] = Int(floor(i/2));
        if (i % 2 == 0)
            d[i] = 1;
        else
            d[i] = 3;
        end
    end
end

function deterministic_form()
    model = Model(GLPK.Optimizer)

    # VARIABLES
    @variable(model, x[n in nodes, t in stage[n]] >= 0, Int)
    @variable(model, w[n in nodes, t in stage[n]] >= 0, Int)
    @variable(model, y[n in nodes, t in stage[n]] >= 0, Int)


    # OBJECTIVE
    @objective(model, Min,
        (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
        (1/2)*sum(x[n,stage[n]] + 3*w[n,stage[n]] + 0.5*y[n,stage[n]] for n=2:3) +
        (1/2)*sum(x[n,stage[n]] + 3*w[n,stage[n]] for n=4:5)
               )

    # CONSTRAINTS
    for n in nodes, t in stage[n]
        # production capacity
        @constraint(model, x[n,t] <= 2)
        # demand
        if stage[n] == 0
            @constraint(model, x[n,t] + w[n,t] - y[n,t] == d[n])
        elseif stage[n] == 1
            @constraint(model, y[1,0] + x[n,t] + w[n,t] - y[n,t] == d[n])
        else
            @constraint(model, y[2,1] + x[n,t] + w[n,t] - y[n,t] == d[n])
            @constraint(model, y[3,1] + x[n,t] + w[n,t] - y[n,t] == d[n])
        end
    end

    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x)
    @show value.(w)
    @show value.(y)
    return
end


function deterministic_form_with_nonant()
    model = Model(GLPK.Optimizer)

    ystages = [[0], [0,1], [0,1,2]];

    # VARIABLES
    @variable(model, x[n in nodes, t in stage[n]] >= 0, Int)
    @variable(model, w[n in nodes, t in stage[n]] >= 0, Int)
    @variable(model, y[n in nodes, t in ystages[stage[n]+1]] >= 0, Int)

    # OBJECTIVE
    @objective(model, Min,
        (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
        (1/2)*sum(x[n,stage[n]] + 3*w[n,stage[n]] + 0.5*y[n,stage[n]] for n=2:3) +
        (1/2)*sum(x[n,stage[n]] + 3*w[n,stage[n]] for n=4:5)
            )

    # CONSTRAINTS
    # non-anticipativity
    for n in 1:(2T-1)-2
        if (n % 2 == 0)
            @constraint(model, [i=2:3], y[n, stage[n]] == y[n+i, stage[n]-1])
        else
            @constraint(model, [i=1:2], y[n, stage[n]] == y[n+i, stage[n]])
        end
    end
    for n in nodes, t in stage[n]
        # production capacity
        @constraint(model, x[n,t] <= 2)
        # demand
        if n == 1
            @constraint(model, x[n,t] + w[n,t] - y[n,t] == d[n])
        else
            @constraint(model, [i=1:stage[n]], y[n,t-i] + x[n,t] + w[n,t] - y[n,t] == d[n])
        end
    end

    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x)
    @show value.(w)
    @show value.(y)
    return
end
