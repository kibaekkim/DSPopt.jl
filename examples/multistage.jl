using MPI
using ScenTrees # This package should be added before the package below. I don't know why...
using StructJuMP
using DSPopt

# scenario tree w/ 3 stages, each parent node has 2 children; dim = 1
tree = Tree([1,2,2],1);
T = 3; # number of stages

# Comment out this line if you want to run in serial
MPI.Init()

# Initialize DSPopt.jl with the communicator.
DSPopt.parallelize(MPI.COMM_WORLD)

d = [[1,0,0,0,0,0,0] [0,1,3,0,0,0,0] [0,0,0,1,3,1,3]] # make this sparse

# create StructuredModel with number of scenarios
model = StructuredModel(num_scenarios = 4)

# VARIABLES
@variable(model, x[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
@variable(model, w[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)
@variable(model, y[n in nodes(tree),t in stage(tree)[n]] >= 0, Int)

@constraint(model, x[1,0] + w[1,0] - y[1,0] == d[1,1]) # demand at root node
@constraint(model, x[1,0] <= 2) # production capacity at root node

@objective(model, Min, 
    x[1,0] + 3*w[1,0] + 0.5*y[1,0]
    + 0.5 * sum(x[k,1] + 3*w[k,1] + 0.5*y[k,1] for k in 2:3) 
    + 0.25 * sum(x[k,2] + 3*w[k,2] for k in 4:7)
)

for m = 2:length(nodes(tree))
    # create a StructuredModel linked to model with id m
    blk = StructuredModel(parent = model, id = m-1)

    # OBJECTIVE
    @objective(blk, Min, 0.)

    # CONSTRAINTS
    for t in stage(tree)[m]
        println("x[$m,$t]")
        println("y[$(root(tree,m)[t]),$(t-1)] + x[$m,$t] + w[$m,$t] - y[$m,$t] == d[$m,$(t+1)]")
        @constraint(blk, x[m,t] <= 2) # production capacity
        @constraint(blk, y[root(tree,m)[t],t-1] + x[m,t] + w[m,t] - y[m,t] == d[m,t+1]) # demand
    end
end

status = optimize!(model,
    is_stochastic = false, # Needs to indicate that the model is a stochastic program.
    solve_type = DSPopt.DW, # see instances(DSPopt.Methods) for other methods
)

if status == MOI.OPTIMAL
    @show objective_value(model)
end

# Comment out this line if you want to run in serial
MPI.Finalize()
