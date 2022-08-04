using ScenTrees
using JuMP
using GLPK

# scenario tree w/ 3 stages, each parent node has 2 children; dim = 1
tree = Tree([1,2,2],1);
T = 3; # number of stages

d = [1, 1, 3, 1, 3, 1, 3]; # demand


function deterministic_form()
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
    model = Model(GLPK.Optimizer)

    ystages = [[0], [0,1], [1,2]];

    # VARIABLES
    @variable(model, x[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
    @variable(model, w[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
    @variable(model, y[n in 1:length(nodes(tree)),t=ystages[stage(tree)[n]+1]] >= 0, Int)

    # OBJECTIVE
    @objective(model, Min,
        (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
        (1/2)*sum(
            x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] + 0.5*y[n,stage(tree)[n]]
                for n=2:3) +
        (1/4)*sum(x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] for n=4:7)
            )

    # CONSTRAINTS
    # non-anticipativity
    @constraint(model, [n=2:length(nodes(tree))],
        y[Int(floor(n/2)), stage(tree)[n]-1] == y[n, stage(tree)[n]-1])
    # production capacity
    @constraint(model, [n=1:length(nodes(tree)), t=stage(tree)[n]], x[n,t] <= 2)
    # demand
    @constraint(model, x[1,0] + w[1,0] - y[1,0] == d[1])
    @constraint(model, [n=2:length(nodes(tree)), t=stage(tree)[n]],
        y[n,t-1] + x[n,t] + w[n,t] - y[n,t] == d[n])

    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x)
    @show value.(w)
    @show value.(y)
    return
end
