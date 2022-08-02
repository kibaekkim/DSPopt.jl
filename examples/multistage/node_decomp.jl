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
        x[1,0] + 3*w[1,0] + 0.5*y[1,0]
        + sum(0.5*(x[k,1] + 3*w[k,1] + 0.5*y[k,1]) for k in 2:3)
        + sum(0.25*(x[k,2] + 3*w[k,2]) for k in 4:7)
    )

    # CONSTRAINTS
    for n = 1:length(nodes(tree))
        # production capacity
        @constraint(model, [t in stage(tree)[n]], x[n,t] <= 2)
        # demand
        if n == 1
            @constraint(model, x[n,0] + w[n,0] - y[n,0] == d[n])
        else
            @constraint(model, [t in stage(tree)[n]],
                y[root(tree,n)[t],t-1] + x[n,t] + w[n,t] - y[n,t] == d[n])
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

    # VARIABLES
    @variable(model, x[n in nodes(tree),t=0:T-1] >= 0, Int)
    @variable(model, w[n in nodes(tree),t=0:T-1] >= 0, Int)
    @variable(model, y[n in nodes(tree),t=0:T-1] >= 0, Int)

    # OBJECTIVE
    @objective(model, Min,
        (1/3) * sum(
            (x[1,t] + 3*w[1,t] + 0.5*y[1,t]) +
            (1/2)*sum(x[n,t] + 3*w[n,t] + 0.5*y[n,t] for n=2:3) +
            (1/4)*sum(x[n,t] + 3*w[n,t] for n=4:7) for t=0:T-1
        )
    )

    # CONSTRAINTS
    for n = 1:length(nodes(tree))
        # production capacity
        @constraint(model, [t=0:T-1], x[n,t] <= 2)
        # non-anticipativity
        @constraint(model, [t=1:T-1], x[n,t-1] == x[n,t])
        @constraint(model, [t=1:T-1], w[n,t-1] == w[n,t])
        @constraint(model, [t=1:T-1], y[n,t-1] == y[n,t])
        # demand
        if n == 1
            @constraint(model, [t=0:T-1], x[n,t] + w[n,t] - y[n,t] == d[n])
        elseif n in 2:3
            @constraint(model, [t=1:T-1], y[root(tree,n)[1],0] + x[n,t] + w[n,t] - y[n,t] == d[n])
        else
            @constraint(model, [t=1:T-1], y[root(tree,n)[2],1] + x[n,t] + w[n,t] - y[n,t] == d[n])
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
