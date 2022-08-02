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

    # OBJECTIVE
    @objective(model, Min,
        x[1,0] + 3*w[1,0] + 0.5*y[1,0]
        + sum(0.5*(x[k,1] + 3*w[k,1] + 0.5*y[k,1]) for k in 2:3)
        + sum(0.25*(x[k,2] + 3*w[k,2]) for k in 4:7)
    )

    # CONSTRAINTS
    @constraint(model, x[1,0] + w[1,0] - y[1,0] == d[1]) # demand at root node
    @constraint(model, x[1,0] <= 2) # production capacity at root node

    for m = 2:length(nodes(tree))
        # create a StructuredModel linked to model with id m-1
        blk = StructuredModel(parent = model, id = m-1)

        # OBJECTIVE
        @objective(blk, Min, 0.)

        # CONSTRAINTS
        @constraint(blk, [t=stage(tree)[m]], x[m,t] <= 2) # production capacity
        @constraint(blk, [t=stage(tree)[m]], y[root(tree,m)[t],t-1] + x[m,t] + w[m,t] - y[m,t] == d[m]) # demand
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

    # OBJECTIVE
    @objective(model, Min,
        (1/3) * sum(
            (x[1,t] + 3*w[1,t] + 0.5*y[1,t]) +
            (1/2)*sum(x[n,t] + 3*w[n,t] + 0.5*y[n,t] for n=2:3) +
            (1/4)*sum(x[n,t] + 3*w[n,t] for n=4:7) for t=0:T-1
        )
    )

    # CONSTRAINTS
    @constraint(model, x[1,0] + w[1,0] - y[1,0] == d[1]) # demand at root node
    @constraint(model, x[1,0] <= 2) # production capacity at root node
    # non-anticipativity
    @constraint(model, [n=nodes(tree),t=1:T-1], x[n,t-1] == x[n,t])
    @constraint(model, [n=nodes(tree),t=1:T-1], w[n,t-1] == w[n,t])
    @constraint(model, [n=nodes(tree),t=1:T-1], y[n,t-1] == y[n,t])

    for m = 2:length(nodes(tree))
        # create a StructuredModel linked to model with id m-1
        blk = StructuredModel(parent = model, id = m-1)

        # OBJECTIVE
        @objective(blk, Min, 0.)

        # CONSTRAINTS
        @constraint(blk, [t=0:T-1], x[m,t] <= 2) # production capacity
        if m in 2:3 # demand
            @constraint(blk, y[root(tree,m)[1],0] + x[m,1] + w[m,1] - y[m,1] == d[m])
        elseif m in 4:7
            @constraint(blk, y[root(tree,m)[2],1] + x[m,2] + w[m,2] - y[m,2] == d[m])
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
