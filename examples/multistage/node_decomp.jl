using ScenTrees
using JuMP
using GLPK

d = [1, 1, 3, 1, 3, 1, 3]; # demand

function deterministic_form()

    # scenario tree w/ 3 stages, each parent node has 2 children; dim = 1
    tree = Tree([1,2,2],1);

    model = Model(GLPK.Optimizer)

    # VARIABLES
    @variable(model, x[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
    @variable(model, w[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
    @variable(model, y[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)


    # OBJECTIVE
    @objective(model, Min,
        x[1,0] + 3*w[1,0] + 0.5*y[1,0] +
        (1/2)*sum(
            x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] + 0.5*y[n,stage(tree)[n]]
            for n in 2:3) +
        (1/4)*sum(x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] for n in 4:7)
            )

    # CONSTRAINTS
    # production capacity
    @constraint(model, [n=1:length(nodes(tree)), t=stage(tree)[n]], x[n,t] <= 2)
    # demand
    @constraint(model, x[1,0] + w[1,0] - y[1,0] == d[1])
    @constraint(model, [n=2:length(nodes(tree)), t=stage(tree)[n]],
        y[root(tree,n)[t],t-1] + x[n,t] + w[n,t] - y[n,t] == d[n])

    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x)
    @show value.(w)
    @show value.(y)
    return
end


function deterministic_form_with_nonant()

    T = 3; # number of stages - MUST BE >= 3 !!

    # scenario tree w/ T stages
    nodes = 1:(2^T)-1; # nodes in tree
    stage = Vector{Int64}(undef, (2^T)-1); # stage of each node in tree
    for i in nodes, t in 1:T
        if i == 1 # stage for root node
            stage[i] = 0;
        elseif i in (2^t):2^(t+1)-1 # stage for all other nodes
            stage[i] = t;
        end
    end

    model = Model(GLPK.Optimizer)

    ystages = [[0], [0,1], [1,2]];

    # VARIABLES
    @variable(model, x[n in nodes, t in stage[n]] >= 0, Int)
    @variable(model, w[n in nodes, t in stage[n]] >= 0, Int)
    @variable(model, y[n in 1:(2^T)-1, t=ystages[stage[n]+1]] >= 0, Int)

    # OBJECTIVE
    @objective(model, Min,
        (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
        (1/2)*sum(
            x[n,stage[n]] + 3*w[n,stage[n]] + 0.5*y[n,stage[n]]
                for n=2:3) +
        (1/4)*sum(x[n,stage[n]] + 3*w[n,stage[n]] for n=4:7)
            )

    # CONSTRAINTS
    # non-anticipativity
    @constraint(model, [n=2:(2^T)-1],
        y[Int(floor(n/2)), stage[n]-1] == y[n, stage[n]-1])
    # production capacity
    @constraint(model, [n=1:(2^T)-1, t=stage[n]], x[n,t] <= 2)
    # demand
    @constraint(model, x[1,0] + w[1,0] - y[1,0] == d[1])
    @constraint(model, [n=2:(2^T)-1, t=stage[n]],
        y[n,t-1] + x[n,t] + w[n,t] - y[n,t] == d[n])

    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x)
    @show value.(w)
    @show value.(y)
    return
end
