using DSPopt

include("farmer_models.jl")

m = farmer_blocks()

"""
NOTE: 
  JuMP v0.21 changes the order of indexing variables.
  This affects some of the test lines below, which were based on the older variable indexing.
"""

@testset "Parent model" begin
    @test length(m.variables) == 27
    @test length(m.constraints) == 6
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(m)
    @test length(start) == length(m.constraints) + 1
    @test start[end] == length(index)
    @test start == [0, 2, 4, 6, 8, 10, 12]
    # @test index == [0, 1, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8] # will fail with JuMP v0.21
    @test value == [1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1., -1.]
    @test rlbd == zeros(6)
    @test rubd == zeros(6)
    # @test obj == [
    #     probability[1]*Cost[1], probability[2]*Cost[1], probability[3]*Cost[1],
    #     probability[1]*Cost[2], probability[2]*Cost[2], probability[3]*Cost[2],
    #     probability[1]*Cost[3], probability[2]*Cost[3], probability[3]*Cost[3],
    #     probability[1]*Purchase[1], probability[2]*Purchase[1], probability[3]*Purchase[1],
    #     probability[1]*Purchase[2], probability[2]*Purchase[2], probability[3]*Purchase[2],
    #     -probability[1]*Sell[1], -probability[2]*Sell[1], -probability[3]*Sell[1],
    #     -probability[1]*Sell[2], -probability[2]*Sell[2], -probability[3]*Sell[2],
    #     -probability[1]*Sell[3], -probability[2]*Sell[3], -probability[3]*Sell[3],
    #     -probability[1]*Sell[4], -probability[2]*Sell[4], -probability[3]*Sell[4]] # will fail with JuMP v0.21
    @test clbd == zeros(27)
    @test cubd == zeros(27) .+ Inf
    @test ctype == "IIIIIIIIICCCCCCCCCCCCCCCCCC"
end

@testset "Child model $i" for (i,subm) in m.children
    @test i > 0
    @test length(subm.variables) == 0
    @test length(subm.constraints) == 5
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(subm)
    @test length(start) == length(subm.constraints) + 1
    @test start[end] == length(index)
    @test start == [0, 3, 6, 9, 12, 13]
    # @test index == [0, 3, 6,  0, 9, 15,  3, 12, 18,  6, 21, 24,  21] .+ (i-1) # will fail with JuMP v0.21
    @test value == [1., 1., 1.,  Yield[i,1], 1., -1.,  Yield[i,2], 1., -1.,  Yield[i,3], -1., -1., 1.]
    @test rlbd == [-Inf; Minreq; -Inf]
    @test rubd == [Budget, Inf, Inf, Inf, 6000]
    @test obj == []
    @test clbd == []
    @test cubd == []
    @test ctype == ""
end

@testset "optimize!: $j" for j in [DSPopt.ExtensiveForm, DSPopt.DW]
    status = DSPopt.optimize!(m, solve_type = j, is_stochastic = false)
    @test status == MOI.OPTIMAL
    @test isapprox(objective_value(m), -108390.)
    @test isapprox(dual_objective_value(m), -108390.)

    primsol = value()
    dualsol = dual()
    for s = 1:3
        @test primsol[s] == []
    end
    @test dualsol == []

    x = m[:x]
    y = m[:y]
    w = m[:w]
    @test isapprox(value(x[1,1]), 170.0)
    @test isapprox(value(x[1,2]), 170.0)
    @test isapprox(value(x[1,3]), 170.0)
    @test isapprox(value(x[2,1]), 80.0)
    @test isapprox(value(x[2,2]), 80.0)
    @test isapprox(value(x[2,3]), 80.0)
    @test isapprox(value(x[3,1]), 250.0)
    @test isapprox(value(x[3,2]), 250.0)
    @test isapprox(value(x[3,3]), 250.0)
    @test isapprox(value(y[1,1]), 0.0)
    @test isapprox(value(y[1,2]), 0.0)
    @test isapprox(value(y[1,3]), 0.0)
    @test isapprox(value(y[2,1]), 0.0)
    @test isapprox(value(y[2,2]), 0.0)
    @test isapprox(value(y[2,3]), 48.0)
    @test isapprox(value(w[1,1]), 310.0)
    @test isapprox(value(w[1,2]), 225.0)
    @test isapprox(value(w[1,3]), 140.0)
    @test isapprox(value(w[2,1]), 48.0)
    @test isapprox(value(w[2,2]), 0.0)
    @test isapprox(value(w[2,3]), 0.0)
    @test isapprox(value(w[3,1]), 6000.0)
    @test isapprox(value(w[3,2]), 5000.0)
    @test isapprox(value(w[3,3]), 4000.0)
    @test isapprox(value(w[4,1]), 0.0)
    @test isapprox(value(w[4,2]), 0.0)
    @test isapprox(value(w[4,3]), 0.0)
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
        @test dsp.solve_type == DSPopt.DW
    end
end