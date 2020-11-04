using StructJuMP
using Ipopt

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

    # farmer2.txt: quadratic constraints 
    @constraint(blk, const_quad, w[1]^2 <= 1600)
    if s == 1
        @constraint(blk, const_quad2, 2 * w[1]^2 + 2 * w[2]^2 <= 400)
    elseif s == 2
        @constraint(blk, const_quad2, 5 * w[1] - w[2]^2 -2 * w[3]^2 >= -100)
    end
end


@testset "optimize!: $j" for j in [DSPopt.Legacy, DSPopt.ExtensiveForm] #instances(DSPopt.Methods)
    
    dsp.is_stochastic = true
    dsp.is_quadratic = true
    
    status = DSPopt.optimize!(m, is_stochastic = true, is_quadratic = true, solve_type = j)
    @test DSPopt.termination_status(m) == MOI.OPTIMAL
    @test isapprox(dual_objective_value(m), -44147, rtol=0.1)

    primsol = value()
    dualsol = dual()
    
    # print("Optimal objective value: ", objective_value(m), "\n")
    # print("Optimal primal solution: \n")
    # print(primsol, "\n")
    
    if dsp.solve_type == DSPopt.Legacy
        for k = 0:3
            @test primsol[k] != []
        end
        @test dualsol != []
    else
        @test isapprox(primsol[0], [96.0, 83.0, 321.0])
        if dsp.solve_type != DSPopt.Benders
            @test isapprox(primsol[1], [0.0, 0.0, 10.60425602775674, 9.35680147726053, 6000.0, 1704.0])
            @test isapprox(primsol[2], [0.0, 0.0, 39.999996661715876, 8.999999638262153, 10.464224243925552, 6409.535775756074])
            @test isapprox(primsol[3], [8.0, 40.8, 0.0, 0.0, 5136.0, 0.0])
        end
        @test dualsol == []
        @test isapprox(value(x[1]), 96.0)
        @test isapprox(value(x[2]), 83.0)
        @test isapprox(value(x[3]), 321.0)
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
        @test dsp.is_quadratic == false
        @test dsp.solve_type == DSPopt.Dual
        @test length(dsp.quadConstrs) == 0
        @test length(dsp.linConstrs) == 0
    end
end

# @testset "writeMps" begin
#     DSPopt.writeMps!(m, "farmer", is_stochastic = true)
#     @test isfile("farmer.mps")
# end