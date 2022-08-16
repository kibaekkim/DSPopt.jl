using StructJuMP
using DSPopt
using MPI

T = 3; # number of stages - MUST BE >= 3 !!

# scenario lattice w/ T stages
nodes = 1:2*T-1; # nodes in lattice
stage = Vector{Int64}(undef, 2*T-1); # stage of each node in lattice
d = Vector{Int64}(undef, 2*T-1); # demand vector
for i in nodes
    if i == 1 # stage/demand for root node
        stage[i] = 0;
        d[i] = 1;
    else # stage/demand for all other nodes
        stage[i] = Int(floor(i/2));
        if (i % 2 == 0)
            d[i] = 1;
        else
            d[i] = 3;
        end
    end
end

MPI.Init()
DSPopt.parallelize(MPI.COMM_WORLD)

# create StructuredModel with number of scenarios
model = StructuredModel(num_scenarios = 2^(T-1))

ystages = Vector{Vector{Int64}}(undef, T); # index for y
ystages[1] = [0];
ystages[2] = [0,1];
for t in 3:T
    ystages[t] = [t-3,t-2,t-1];
end

# VARIABLES
@variable(model, x[n in nodes,  t in stage[n]] >= 0, Int)
@variable(model, w[n in nodes,  t in stage[n]] >= 0, Int)
@variable(model, y[n in 1:2T-1, t in ystages[stage[n]+1]] >= 0, Int)

# OBJECTIVE (parent)
@objective(model, Min,
    (x[1,0] + 3*w[1,0] + 0.5*y[1,0]) +
    (1/2)*sum(x[n,stage[n]] + 3*w[n,stage[n]] + 0.5*y[n,stage[n]] for n=2:2*T-3) +
    (1/2)*sum(x[n,stage[n]] + 3*w[n,stage[n]] for n=2*T-2:2*T-1)
           )

# CONSTRAINTS (parent): non-anticipativity
for n in 1:(2*T-1)-2
    if (n % 2 == 0)
        @constraint(model, [i=2:3], y[n, stage[n]] == y[n+i, stage[n]-1])
    else
        @constraint(model, [i=1:2], y[n, stage[n]] == y[n+i, stage[n]])
    end
end

for m = 1:2*T-1
    # create a StructuredModel linked to model with id m
    blk = StructuredModel(parent = model, id = m)

    # OBJECTIVE (subprob)
    @objective(blk, Min, 0.)

    # CONSTRAINTS (subprob)
    # production capacity
    @constraint(blk, [t=stage[m]], x[m,t] <= 2)
    @constraint(blk, [t=stage[m]], w[m,t] <= maximum(d)) # need upper bounds on all vars for some reason
    @constraint(blk, [t=stage[m]], y[m,t] <= maximum(d))

    # demand
    if m == 1
        @constraint(blk, x[m,0] + w[m,0] - y[m,0] == d[m])
    else
        @constraint(blk, [t=stage[m], i=1:stage[m]], y[m,t-i] + x[m,t] + w[m,t] - y[m,t] == d[m])
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
