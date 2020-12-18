using DSPopt

include("farmer_models.jl")

m = farmer_stochastic()

@testset "Parent model" begin
    @test length(m.variables) == 3
    @test length(m.constraints) == 1
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(m)
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
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(subm)
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

@testset "optimize!: $j" for j in instances(DSPopt.Methods)
    dsp.is_stochastic = true
    status = DSPopt.optimize!(m, is_stochastic = true, solve_type = j)
    @test DSPopt.termination_status(m) == MOI.OPTIMAL
    if dsp.solve_type in [DSPopt.DW, DSPopt.ExtensiveForm]
        @test isapprox(objective_value(m), -108390.)
    else
        @test objective_value(m) >= -108390.
    end
    @test isapprox(dual_objective_value(m), -108390., atol=0.1)

    x = m[:x]
    primsol = value()
    dualsol = dual()
    if dsp.solve_type == DSPopt.Dual
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
