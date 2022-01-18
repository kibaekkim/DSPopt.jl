"""
Kibaek Kim - ANL MCS 2016
  Updated in 2020
Farmer example from Birge and Louveaux book.

To use MPI, add 
- MPI.Init()
- DSPopt.parallelize(comm)
- MPI.Finalize()
"""

using MPI
using StructJuMP
using DSPopt

MPI.Init()

# Initialize DSPopt.jl with the communicator.
DSPopt.parallelize(MPI.COMM_WORLD)

NS = 3;                        # number of scenarios
probability = [1/3, 1/3, 1/3]; # probability

CROPS  = 1:3 # set of crops (wheat, corn and sugar beets, resp.)
PURCH  = 1:2 # set of crops to purchase (wheat and corn, resp.)
SELL   = 1:4 # set of crops to sell (wheat, corn, sugar beets under 6K and those over 6K)

Cost     = [150 230 260]   # cost of planting crops
Budget   = 500             # budget capacity
Purchase = [238 210];      # purchase price
Sell     = [170 150 36 10] # selling price
Yield    = [3.0 3.6 24.0;
            2.5 3.0 20.0;
            2.0 2.4 16.0]
Minreq   = [200 240 0]     # minimum crop requirement

m = StructuredModel(num_scenarios = NS)

@variable(m, 0 <= x[i=CROPS] <= 500, Int)
@objective(m, Min, sum(Cost[i] * x[i] for i=CROPS))
@constraint(m, const_budget, sum(x[i] for i=CROPS) <= Budget)

# TODO: Distributing blocks is not supported yet. 
# This is because of the DSPopt implementation.
for s in 1:NS
    blk = StructuredModel(parent = m, id = s, prob = probability[s])

    @variable(blk, y[j=PURCH] >= 0)
    @variable(blk, w[k=SELL] >= 0)

    @objective(blk, Min, sum(Purchase[j] * y[j] for j=PURCH) - sum(Sell[k] * w[k] for k=SELL))

    @constraint(blk, const_minreq[j=PURCH], Yield[s,j] * x[j] + y[j] - w[j] >= Minreq[j])
    @constraint(blk, const_minreq_beets, Yield[s,3] * x[3] - w[3] - w[4] >= Minreq[3])
    @constraint(blk, const_aux, w[3] <= 6000)
end

status = optimize!(m, 
    is_stochastic = true, # Needs to indicate that the model is of the stochastic program.
    solve_type = DSPopt.DW, # see instances(DSPopt.Methods) for other methods
    param = "examples/params.txt" # This path assumes running from the one-level upper directory (i.e., ../).
    )

if DSPopt.myrank() == 0 && status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
    @show objective_value(m)
    @show dual_objective_value(m)
    @show value.(x)
    @show dual() # This is available only for solve_type = DSPopt.Dual.
end

MPI.Finalize()
