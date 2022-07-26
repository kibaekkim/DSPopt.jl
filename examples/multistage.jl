using Pkg
Pkg.add("ScenTrees")
Pkg.add("StructJuMP")
Pkg.add("DSPopt")
Pkg.activate(".");

using ScenTrees # This package should be added before the package below. I don't know why...
using StructJuMP
using DSPopt

# scenario tree w/ 3 stages, each parent node has 2 children; dim = 1
tree = Tree([1,2,2],1);
T = 3; # number of stages

d = [1,1,3,1,3,1,3] # demand vector

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
# assuming overtime labor is a corrective measure that can be taken once uncertainty is realized
for r in 2:length(nodes(tree))-1, s in r+1:length(nodes(tree)) # non-anticipativity
    if root(tree,r)[stage(tree)[r]] == root(tree,s)[stage(tree)[s]]
        @constraint(model, x[r,stage(tree)[r]] == x[s,stage(tree)[s]])
        #@constraint(model, w[r,stage(tree)[r]] == w[s,stage(tree)[s]])
    end
end

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
