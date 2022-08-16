using JuMP
using GLPK

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

d = Vector{Int64}(undef, (2^T)-1); # demand vector
for i in nodes
    if i == 1 # demand for root node
        d[i] = 1;
    else # demand for all other nodes
        if (i % 2 == 0)
            d[i] = 1;
        else
            d[i] = 3;
        end
    end
end

model = Model(GLPK.Optimizer)

ystages = Vector{Vector{Int64}}(undef, T); # index for y
ystages[1] = [0];
for t in 2:T
    ystages[t] = [t-2,t-1];
end

# VARIABLES
@variable(model, x[n in nodes, t in stage[n]] >= 0, Int)
@variable(model, w[n in nodes, t in stage[n]] >= 0, Int)
@variable(model, y[n in 1:(2^T)-1, t=ystages[stage[n]+1]] >= 0, Int)

# OBJECTIVE
@objective(model, Min,
    (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
    sum((1/(2^(stage[n])))*(x[n,stage[n]] + 3*w[n,stage[n]] + 0.5*y[n,stage[n]]) for n in 2:2^(T-1)-1) +
    sum((1/(2^(stage[n])))*(x[n,stage[n]] + 3*w[n,stage[n]]) for n in 2^(T-1):(2^T)-1)
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
