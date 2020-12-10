using DSPopt
using Test
using SparseArrays

const dsp = DSPopt.dspenv

@testset "Initializing DSPopt" begin
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
    @test isnothing(dsp.comm)
    @test dsp.comm_size == 1
    @test dsp.comm_rank == 0
end

@testset "Setting options" begin
    @testset "param:" begin
        DSPopt.setoptions!(Dict(:param => "params.txt"))
    end
    @testset "is_stochastic: $i" for i in [true, false]
        DSPopt.setoptions!(Dict(:is_stochastic => i))
        @test dsp.is_stochastic == i
    end
    @testset "solve_type: $t" for t in instances(DSPopt.Methods)
        DSPopt.setoptions!(Dict(:solve_type => t))
        @test dsp.solve_type == t
    end
end

@testset "Farmer example: stochastic form" begin
    include("farmer_stoc.jl")
end

@testset "Farmer example: block form" begin
    include("farmer_block.jl")
end

@testset "dcap" begin
    include("dcap.jl")
end

@testset "Farmer example: stochastic quadratic form" begin
    include("farmer_qcp.jl")
end

@testset "Farmer example2: stochastic quadratic form" begin
    include("farmer_qcp2.jl")
end

@testset "Farmer example3: stochastic quadratic form" begin
    include("farmer_qcp3.jl")
end

@testset "Freeing DSPopt" begin
    DSPopt.freeEnv(dsp)
    @test dsp.p == C_NULL
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
    @test isnothing(dsp.comm)
    @test dsp.comm_size == 1
    @test dsp.comm_rank == 0
end