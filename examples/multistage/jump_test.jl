using JuMP
using GLPK

d = [1,1,3,1,3,1,3] # demand vector

model = Model(GLPK.Optimizer)

# VARIABLES
@variable(model, x[1:7] >= 0, Int)
@variable(model, w[1:7] >= 0, Int)
@variable(model, y[1:7] >= 0, Int)

# OBJECTIVE
@objective(model, Min,
    x[1] + 3*w[1] + 0.5*y[1]
    + sum(0.5*(x[k] + 3*w[k] + 0.5*y[k]) for k in 2:3)
    + sum(0.25*(x[k] + 3*w[k]) for k in 4:7)
)

# CONSTRAINTS
@constraint(model, x[1] + w[1] - y[1] == d[1]) # demand at root node
@constraint(model, [n=1:7], x[n] <= 2) # production capacity
# non-anticipativity
@constraint(model, x[2] == x[3])
@constraint(model, x[4] == x[5])
@constraint(model, x[6] == x[7])
#demand
@constraint(model, [n=2:3], y[1] + x[n] + w[n] - y[n] == d[n])
@constraint(model, [n=4:5], y[2] + x[n] + w[n] - y[n] == d[n])
@constraint(model, [n=6:7], y[3] + x[n] + w[n] - y[n] == d[n])

optimize!(model)
termination_status(model)
value.(x)
value.(w)
value.(y)
