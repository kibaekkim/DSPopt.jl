# DSPopt.jl
![Run tests](https://github.com/kibaekkim/DSPopt.jl/workflows/Run%20tests/badge.svg)
[![codecov](https://codecov.io/gh/kibaekkim/DSPopt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kibaekkim/DSPopt.jl)

DSPopt.jl is an interface to a parallel decomposition mixed-integer programming solver [DSP](https://github.com/Argonne-National-Laboratory/DSP). 
This package allows users to define block structures in optimization model written in [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl) 
and solve the block-structured problem using the parallle solver ``DSP``.

This package can model and solve (by interfacing `DSP`) the following optimization problems:

- Two-stage stochastic (mixed-integer) quadratic constrained programming (TSSP)
- The Wasserstein-based distributionally robust variant of TSSP
- Block-structured (mixed-integer) linear programming

## Intallation

> **_NOTE:_** You need to install solver [DSP](https://github.com/Argonne-National-Laboratory/DSP) first. This package provides an interface only.

```julia
] add DSPopt
```

## Examples

This is one simple example of `stochastic` form.
Please find more examples in `./examples` particularly for `block` form.

```julia
using MPI
using StructJuMP
using DSPopt

# Comment out this line if you want to run in serial
MPI.Init()

# Initialize DSPopt.jl with the communicator.
DSPopt.parallelize(MPI.COMM_WORLD)

xi = [[7,7] [11,11] [13,13]]

# create StructuredModel with number of scenarios
m = StructuredModel(num_scenarios = 3)

@variable(m, 0 <= x[i=1:2] <= 5, Int)
@objective(m, Min, -1.5 * x[1] - 4 * x[2])
for s = 1:3
    # create a StructuredModel linked to m with id s and probability 1/3
    blk = StructuredModel(parent = m, id = s, prob = 1/3)
    @variable(blk, y[j=1:4], Bin)
    @objective(blk, Min, -16 * y[1] + 19 * y[2] + 23 * y[3] + 28 * y[4])
    @constraint(blk, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= xi[1,s] - x[1])
    @constraint(blk, 6 * y[1] + y[2] + 3 * y[3] + 2 * y[4] <= xi[2,s] - x[2])

    # Quadratic constraints supported from DSP version 2.0.0 or higher
    @constraint(blk, const_quad, w[1]^2 <= 1600)
    if s == 1
        @constraint(blk, const_quad2, 2 * w[1]^2 + 2 * w[2]^2 <= 400)
    elseif s == 2
        @constraint(blk, const_quad2, 5 * w[1] - w[2]^2 -2 * w[3]^2 >= -100)
    end
end

# The Wasserstein ambiguity set of order-2 can be imposed with the size limit of 1.0.
DSPopt.set(WassersteinSet, 2, 1.0)

status = optimize!(m, 
    is_stochastic = true, # Needs to indicate that the model is a stochastic program.
    solve_type = DSPopt.DW, # see instances(DSPopt.Methods) for other methods
)

# Comment out this line if you want to run in serial
MPI.Finalize()
```

## Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.
