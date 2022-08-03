using ScenTrees # This package should be added before the package below. I don't know why...
using StructJuMP
using DSPopt

# scenario tree w/ 3 stages, each parent node has 2 children; dim = 1
tree = Tree([1,2,2],1);
T = 3; # number of stages

d = [1,1,3,1,3,1,3] # demand vector


function deterministic_form()

    # create StructuredModel with number of scenarios
    model = StructuredModel(num_scenarios = length(nodes(tree))-1)

    # VARIABLES
    @variable(model, x[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
    @variable(model, w[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
    @variable(model, y[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)

    # OBJECTIVE (parent)
    @objective(model, Min,
        (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
        (1/2)*sum(
            x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] + 0.5*y[n,stage(tree)[n]]
            for n in 2:3) +
        (1/4)*sum(x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] for n in 4:7)
            )

    for m = 1:length(nodes(tree))
        # create a StructuredModel linked to model with id m
        blk = StructuredModel(parent = model, id = m)

        # OBJECTIVE (subprob)
        @objective(blk, Min, 0.)

        # CONSTRAINTS (subprob)
        # production capacity
        @constraint(blk, [t=stage(tree)[m]], x[m,t] <= 2)
        # demand
        if m == 1
            @constraint(model, x[m,0] + w[m,0] - y[m,0] == d[m])
        else
            @constraint(blk, [t=stage(tree)[m]],
                y[root(tree,m)[t],t-1] + x[m,t] + w[m,t] - y[m,t] == d[m])
        end
    end

    status = optimize!(model,
        is_stochastic = false, # Needs to indicate that the model is NOT a stochastic program.
        solve_type = DSPopt.ExtensiveForm, # see instances(DSPopt.Methods) for other methods
    )

    if status == MOI.OPTIMAL
        @show objective_value(model)
        @show value.(x)
        @show value.(w)
        @show value.(y)
    end
    return
end


function deterministic_form_with_nonant()

    # create StructuredModel with number of scenarios
    model = StructuredModel(num_scenarios = length(nodes(tree))-1)

    # VARIABLES
    @variable(model, x[n in nodes(tree),t=0:T-1] >= 0, Int)
    @variable(model, w[n in nodes(tree),t=0:T-1] >= 0, Int)
    @variable(model, y[n in nodes(tree),t=0:T-1] >= 0, Int)

    # OBJECTIVE (parent)
    @objective(model, Min,
        (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
        (1/2)*sum(
            x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] + 0.5*y[n,stage(tree)[n]]
                for n=2:3) +
        (1/4)*sum(x[n,stage(tree)[n]] + 3*w[n,stage(tree)[n]] for n=4:7)
            )

    # CONSTRAINTS (parent)
    # non-anticipativity
    @constraint(model, [n=2:length(nodes(tree))],
        y[Int(floor(n/2)), stage(tree)[n]-1] == y[n, stage(tree)[n]-1])

    for m = 1:length(nodes(tree))
        # create a StructuredModel linked to model with id m
        blk = StructuredModel(parent = model, id = m)

        # OBJECTIVE (subprob)
        @objective(blk, Min, 0.)

        # CONSTRAINTS (subprob)
        # production capacity
        @constraint(blk, [t=0:T-1], x[m,t] <= 2)
        # demand
        if m == 1
            @constraint(blk, x[m,0] + w[m,0] - y[m,0] == d[m])
        else
            @constraint(blk, [t=stage(tree)[m]],
                y[m,t-1] + x[m,t] + w[m,t] - y[m,t] == d[m])
        end
    end

    status = optimize!(model,
        is_stochastic = false, # Needs to indicate that the model is NOT a stochastic program.
        solve_type = DSPopt.ExtensiveForm, # see instances(DSPopt.Methods) for other methods
    )

    if status == MOI.OPTIMAL
        @show objective_value(model)
        @show value.(x)
        @show value.(w)
        @show value.(y)
    end
    return
end
