using DSPopt
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

# JuMP model
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

@testset "Parent model" begin
    @test length(m.variables) == 3
    @test length(m.constraints) == 1
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(m, m.constraints)
    @test length(start) == length(m.constraints) + 1
    @test start[end] == length(index)
    @test start == [0, 3]
    @test index == [0, 1, 2]
    @test value == [1., 1., 1.]
    @test rlbd == [-Inf]
    @test rubd == [Budget]
    @test obj == Cost
    @test clbd == zeros(3)
    @test cubd == zeros(3) .+ Inf
    @test ctype == "III"
end

@testset "Child model $i" for (i,subm) in m.children
    @test i > 0
    @test m.probability[i] == probability[i]
    @test length(subm.variables) == 6
    @test length(subm.constraints) == 4
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(subm, subm.constraints)
    @test length(start) == length(subm.constraints) + 1
    @test start[end] == length(index)
    @test start == [0, 3, 6, 9, 10]
    @test index == [0, 3, 5, 1, 4, 6, 2, 7, 8, 7]
    @test value == [Yield[i,1], 1., -1., Yield[i,2], 1., -1., Yield[i,3], -1., -1., 1.]
    @test rlbd == [Minreq; -Inf]
    @test rubd == [Inf, Inf, Inf, 6000]
    @test obj == [Purchase; -Sell]
    @test clbd == zeros(6)
    @test cubd == zeros(6) .+ Inf
    @test ctype == "CCCCCC"
end

@testset "optimize!: $j" for j in [DSPopt.Legacy, DSPopt.ExtensiveForm]#instances(DSPopt.Methods)
    dsp.is_stochastic = true
    status = DSPopt.optimize!(m, is_stochastic = true, solve_type = j)
    @test DSPopt.termination_status(m) == MOI.OPTIMAL
    if dsp.solve_type in [DSPopt.Dual, DSPopt.ExtensiveForm]
        @test isapprox(objective_value(m), -108390.)
    else
        @test objective_value(m) >= -108390.
    end
    @test isapprox(dual_objective_value(m), -108390., atol=0.1)

    primsol = value()
    dualsol = dual()
    if dsp.solve_type == DSPopt.Legacy
        for k = 0:3
            @test primsol[k] != []
        end
        @test dualsol != []
    else
        @test isapprox(primsol[0], [170.0, 80.0, 250.0])
        if dsp.solve_type != DSPopt.Benders
            @test isapprox(primsol[1], [0.0, 0.0, 310.0, 48.0, 6000.0, 0.0])
            @test isapprox(primsol[2], [0.0, 0.0, 225.0, 0.0, 5000.0, 0.0])
            @test isapprox(primsol[3], [0.0, 48.0, 140.0, 0.0, 4000.0, 0.0])
        end
        @test dualsol == []
        @test isapprox(value(x[1]), 170.0)
        @test isapprox(value(x[2]), 80.0)
        @test isapprox(value(x[3]), 250.0)
    end
    DSPopt.freeSolver(dsp)
    @testset "freeModel" begin
        DSPopt.freeModel(dsp)
        @test dsp.p != C_NULL
        @test length(dsp.numRows) == 0
        @test length(dsp.numCols) == 0
        @test isnan(dsp.primVal)
        @test isnan(dsp.dualVal)
        @test length(dsp.colVal) == 0
        @test length(dsp.rowVal) == 0
        @test dsp.nblocks == -1
        @test dsp.block_ids == []
        @test dsp.is_stochastic == false
        @test dsp.solve_type == DSPopt.Dual
        @test length(dsp.quadConstrs) == 0
        @test length(dsp.linConstrs) == 0
    end
end

@testset "writeMps" begin
    DSPopt.writeMps!(m, "farmer", is_stochastic = true)
    @test isfile("farmer.mps")
end