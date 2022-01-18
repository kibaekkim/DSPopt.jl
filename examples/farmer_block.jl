"""
Kibaek Kim - ANL MCS 2016
  Updated in 2020
Farmer example from Birge and Louveaux book.

This example writes farmer problem in a non-stochastic form 
by explictly defining nonanticipativity constraints on x.
This form particularly represents a block-angular structure of the model.
In this form, each scenario index is interpreted to represent a sub-block.
All variables should be declared at the parent model.
"""

using StructJuMP
using DSPopt

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

@variable(m, x[i=CROPS,s=1:NS] >= 0, Int)
@variable(m, y[j=PURCH,s=1:NS] >= 0)
@variable(m, w[k=SELL,s=1:NS] >= 0)
@objective(m, Min, 
	  sum(probability[s] * Cost[i] * x[i,s] for i=CROPS for s=1:NS)
    + sum(probability[s] * Purchase[j] * y[j,s] for j=PURCH for s=1:NS) 
	- sum(probability[s] * Sell[k] * w[k,s] for k=SELL for s=1:NS))
@constraint(m, nonant[i=CROPS,s=2:NS], x[i,s-1] - x[i,s] == 0)

for s in 1:NS
    blk = StructuredModel(parent = m, id = s)
	@objective(blk, Min, 0.)
	@constraint(blk, const_budget, sum(x[i,s] for i=CROPS) <= Budget)
    @constraint(blk, const_minreq[j=PURCH], Yield[s,j] * x[j,s] + y[j,s] - w[j,s] >= Minreq[j])
    @constraint(blk, const_minreq_beets, Yield[s,3] * x[3,s] - w[3,s] - w[4,s] >= Minreq[3])
    @constraint(blk, const_aux, w[3,s] <= 6000)
end

status = optimize!(m, 
    is_stochastic = false, # Needs to indicate that the model is NOT a stochastic program.
    solve_type = DSPopt.ExtensiveForm, # see instances(DSPopt.Methods) for other methods
    param = "examples/params.txt" # This path assumes running from the one-level upper directory (i.e., ../).
    )

if status in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
    @show objective_value(m)
    @show dual_objective_value(m)
    @show value.(x)
    @show dual() # This is available only for solve_type = DSPopt.Dual.
end
