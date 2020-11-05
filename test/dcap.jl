#=
Source:
  S. Ahmed and R. Garcia. "Dynamic Capacity Acquisition and Assignment under Uncertainty," Annals of Operations Research, vol.124, pp. 267-283, 2003

Input:
  nR: number of resources
  nN: number of tasks
  nT: number of time periods
  nS: number of scenarios

Sets:
  R: resources
  N: tasks
  T: time periods
  S: scenarios

Variables (1st Stage):
  x[i,t]: capacity acquired for resource i at period t
  u[i,t]: 1 if x[i,t] > 0, 0 otherwise

Variables (2nd Stage):
  y[i,j,t]: 1 if resource i is assigned to task j in period t, 0 otherwise

Parameters (general):
  a[i,t]: linear component of expansion cost for resource i at period t
  b[i,t]: fixed component of expansion cost for resource i at period t
  c[i,j,t,s]: cost of assigning resource i to task j in period t
  c0[j,t,s]: penalty incurred if task j in period t is not served

Parameters (scenario):
  d[j,t,s]: capacity required for to perform task j in period t in scenario s
=#

using StructJuMP
using Random

function DCAP(nR::Int, nN::Int, nT::Int, nS::Int, seed::Int=1)::StructuredModel

    # set random seed (default=1)
    Random.seed!(seed)

    # generate & store instance data
    ## sets
    R = 1:nR
    N = 1:nN
    T = 1:nT
    S = 1:nS

    ## parameters
    a = rand(nR, nT) * 5 .+ 5
    b = rand(nR, nT) * 40 .+ 10
    c = rand(nR, nN, nT, nS) * 5 .+ 5
    c0 = rand(nN, nT, nS) * 500 .+ 500
    d = rand(nN, nT, nS) .+ 0.5
    Pr = ones(nS)/nS

    # construct JuMP.Model
    model = StructuredModel(num_scenarios = nS)

    ## 1st stage
    @variable(model, x[i=R,t=T] >= 0)
    @variable(model, u[i=R,t=T], Bin)
    @objective(model, Min, sum(a[i,t]*x[i,t] + b[i,t]*u[i,t] for i in R for t in T))
    @constraint(model, [i=R,t=T], x[i,t] - u[i,t] <= 0)

    ## 2nd stage
    for s in S
        sb = StructuredModel(parent=model, id = s, prob = Pr[s])
        @variable(sb, y[i=R, j=N, t=T], Bin)
        #@variable(sb, z[j=N,t=T] >= 0) # originally implemented variable (continuous)
        @variable(sb, z[j=N,t=T], Bin)  # modify as SIPLIB 1.0
        @objective(sb, Min, sum(c[i,j,t,s]*y[i,j,t] for i in R for j in N for t in T) + sum(c0[j,t,s]*z[j,t] for j in N for t in T))
        @constraint(sb, [i=R, t=T], -sum(x[i,tau] for tau in 1:t) + sum(d[j,t,s]*y[i,j,t] for j in N) <= 0)
        @constraint(sb, [j=N, t=T], sum(y[i,j,t] for i in R) + z[j,t] == 1)
    end

    return model
end

m = DCAP(2,2,3,3);

@testset "Parent model" begin
  @test length(m.variables) == 12
  @test length(m.constraints) == 6
  start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(m, m.constraints)
  @show clbd
  @show cubd
  @show ctype
  # print(m)
end

@testset "Child model $i" for (i,subm) in m.children
  start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSPopt.get_model_data(subm, subm.constraints)
  @show clbd
  @show cubd
  @show ctype
end

@testset "write MPS" begin
  DSPopt.writeMps!(m, "dcap")
  @test isfile("dcap.mps")
end

@testset "optimize!: $j" for j in [DSPopt.ExtensiveForm]
  status = DSPopt.optimize!(m, solve_type = j, is_stochastic = true, param = "params.txt")
end
