using StructJuMP

include("farmer_models.jl")

m = farmer_stochastic()
DSPopt.classifyConstrs(m)

@testset "Parent model" begin
    @test length(m.variables) == 3
    @test length(m.constraints) == 1
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval = DSPopt.get_model_data(m)
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
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname, nqrows, linnzcnt, quadnzcnt, rhs, sense, linstart, linind, linval, quadstart, quadrow, quadcol, quadval = DSPopt.get_model_data(subm, i)
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

DSPopt.set(WassersteinSet, 2, 1.0)

@testset "optimize!: $j" for j in [DSPopt.Dual, DSPopt.Benders]
    status = DSPopt.optimize!(m, is_stochastic = true, solve_type = j)
    @test DSPopt.termination_status(m) == MOI.OPTIMAL
    DSPopt.freeSolver(dsp)
    DSPopt.freeModel(dsp)
end
