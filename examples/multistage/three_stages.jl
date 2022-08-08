using ScenTrees # This package should be added before the package below. I don't know why...
using StructJuMP
using DSPopt
using MPI

# scenario tree w/ 3 stages, each parent node has 2 children; dim = 1
tree = Tree([1,2,2],1);
T = 3; # number of stages

d = [1,1,3,1,3,1,3] # demand vector

MPI.Init()
DSPopt.parallelize(MPI.COMM_WORLD)

# create StructuredModel with number of scenarios
model = StructuredModel(num_scenarios = length(leaves(tree)[2]))

ystages = [[0], [0,1], [1,2]];

# VARIABLES
@variable(model, x[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
@variable(model, w[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
@variable(model, y[n in 1:length(nodes(tree)),t=ystages[stage(tree)[n]+1]] >= 0, Int)

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
    @constraint(blk, [t=stage(tree)[m]], x[m,t] <= 2)
    @constraint(blk, [t=stage(tree)[m]], w[m,t] <= maximum(d)) # need upper bounds on all vars for some reason
    @constraint(blk, [t=stage(tree)[m]], y[m,t] <= maximum(d))

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
    solve_type = DSPopt.DW, # see instances(DSPopt.Methods) for other methods
)

if DSPopt.myrank() == 0 && status == MOI.OPTIMAL
    @show objective_value(model)
    @show value.(x)
    @show value.(w)
    @show value.(y)
end

MPI.Finalize()
