using StructJuMP

NS = 3;                        # number of scenarios
probability = [1/3, 1/3, 1/3]; # probability

CROPS  = 1:3 # set of crops (wheat, corn and sugar beets, resp.)
PURCH  = 1:2 # set of crops to purchase (wheat and corn, resp.)
SELL   = 1:4 # set of crops to sell (wheat, corn, sugar beets under 6K and those over 6K)

Cost     = [150., 230., 260.]    # cost of planting crops
Budget   = 500.                # budget capacity
Purchase = [238., 210.];        # purchase price
Sell     = [170., 150., 36., 10.] # selling price
Yield    = [3.0 3.6 24.0;
            2.5 3.0 20.0;
            2.0 2.4 16.0]
Minreq   = [200., 240., 0.]      # minimum crop requirement

function farmer_stochastic()
    m = StructuredModel(num_scenarios = NS)

    @variable(m, x[i=CROPS] >= 0, Int)
    @objective(m, Min, sum(Cost[i] * x[i] for i=CROPS))
    @constraint(m, sum(x[i] for i=CROPS) <= Budget)

    for s in 1:NS
        blk = StructuredModel(parent = m, id = s, prob = probability[s])
        @variable(blk, y[j=PURCH] >= 0)
        @variable(blk, w[k=SELL] >= 0)
        @objective(blk, Min, 
            sum(Purchase[j] * y[j] for j=PURCH) 
            - sum(Sell[k] * w[k] for k=SELL))
        @constraint(blk, const_minreq[j=PURCH], Yield[s,j] * x[j] + y[j] - w[j] >= Minreq[j])
        @constraint(blk, const_minreq_beets, Yield[s,3] * x[3] - w[3] - w[4] >= Minreq[3])
        @constraint(blk, const_aux, w[3] <= 6000)
    end

    return m
end

function farmer_blocks()
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
        @constraint(blk, sum(x[i,s] for i=CROPS) <= Budget)
        @constraint(blk, const_minreq[j=PURCH], Yield[s,j] * x[j,s] + y[j,s] - w[j,s] >= Minreq[j])
        @constraint(blk, const_minreq_beets, Yield[s,3] * x[3,s] - w[3,s] - w[4,s] >= Minreq[3])
        @constraint(blk, const_aux, w[3,s] <= 6000)
    end
    return m
end