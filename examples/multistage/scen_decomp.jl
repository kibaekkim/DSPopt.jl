using JuMP
using GLPK

d_1 = 1; # demand in 1st stage
d_2 = [1,3]; # demand vector for 2nd stage vars
d_3 = [1,3,1,3]; # demand vector for 3rd stage vars

function deterministic_form()
    model = Model(GLPK.Optimizer)

    # VARIABLES
    @variable(model, x_1 >= 0, Int)
    @variable(model, w_1 >= 0, Int)
    @variable(model, y_1 >= 0, Int)
    @variable(model, x_2[k=1:2] >= 0, Int)
    @variable(model, w_2[k=1:2] >= 0, Int)
    @variable(model, y_2[k=1:2] >= 0, Int)
    @variable(model, x_3[k=1:4] >= 0, Int)
    @variable(model, w_3[k=1:4] >= 0, Int)
    @variable(model, y_3[k=1:4] >= 0, Int)

    # OBJECTIVE
    @objective(model, Min,
        x_1 + 3*w_1 + 0.5*y_1
        + sum(0.5*(x_2[k] + 3*w_2[k] + 0.5*y_2[k]) for k in 1:2)
        + sum(0.25*(x_3[k] + 3*w_3[k]) for k in 1:4)
    )

    # CONSTRAINTS
    @constraint(model, x_1 + w_1 - y_1 == d_1) # demand at root node
    @constraint(model, x_1 <= 2) # production capacity
    @constraint(model, [k=1:2], x_2[k] <= 2) # production capacity
    @constraint(model, [k=1:4], x_3[k] <= 2) # production capacity

    #demand
    @constraint(model, [k=1:2], y_1 + x_2[k] + w_2[k] - y_2[k] == d_2[k])
    @constraint(model, [k=1:2], y_2[1] + x_3[k] + w_3[k] - y_3[k] == d_3[k])
    @constraint(model, [k=3:4], y_2[2] + x_3[k] + w_3[k] - y_3[k] == d_3[k])

    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x_1)
    @show value.(w_1)
    @show value.(y_1)
    @show value.(x_2)
    @show value.(w_2)
    @show value.(y_2)
    @show value.(x_3)
    @show value.(w_3)
    @show value.(y_3)
    return
end

function deterministic_form_with_nonant()
    model = Model(GLPK.Optimizer)

    # VARIABLES
    @variable(model, x_1[k=1:4] >= 0, Int)
    @variable(model, w_1[k=1:4] >= 0, Int)
    @variable(model, y_1[k=1:4] >= 0, Int)
    @variable(model, x_2[k=1:4] >= 0, Int)
    @variable(model, w_2[k=1:4] >= 0, Int)
    @variable(model, y_2[k=1:4] >= 0, Int)
    @variable(model, x_3[k=1:4] >= 0, Int)
    @variable(model, w_3[k=1:4] >= 0, Int)
    @variable(model, y_3[k=1:4] >= 0, Int)

    # OBJECTIVE
    @objective(model, Min,
        0.25 * sum(
            x_1[k] + 3*w_1[k] + 0.5*y_1[k] 
            + (x_2[k] + 3*w_2[k] + 0.5*y_2[k])
            + (x_3[k] + 3*w_3[k])
            for k=1:4
        )
    )

    # CONSTRAINTS
    for k = 1:4
        @constraint(model, x_1[k] + w_1[k] - y_1[k] == d_1) # demand at root node
        @constraint(model, x_1[k] <= 2) # production capacity
        @constraint(model, x_2[k] <= 2) # production capacity
        @constraint(model, x_3[k] <= 2) # production capacity

        #demand
        @constraint(model, y_1[k] + x_2[k] + w_2[k] - y_2[k] == d_2[Int(ceil(k/2))])
        @constraint(model, y_2[k] + x_3[k] + w_3[k] - y_3[k] == d_3[k])
        @constraint(model, y_2[k] + x_3[k] + w_3[k] - y_3[k] == d_3[k])
    end

    # non-anticipativity
    @constraint(model, [k=2:4], x_1[k-1] == x_1[k])
    @constraint(model, [k=2:4], w_1[k-1] == w_1[k])
    @constraint(model, [k=2:4], y_1[k-1] == y_1[k])
    @constraint(model, [k=[1,3]], x_2[k] == x_2[k+1])
    @constraint(model, [k=[1,3]], w_2[k] == w_2[k+1])
    @constraint(model, [k=[1,3]], y_2[k] == y_2[k+1])

    optimize!(model)
    @show termination_status(model)
    @show objective_value(model)
    @show value.(x_1)
    @show value.(w_1)
    @show value.(y_1)
    @show value.(x_2)
    @show value.(w_2)
    @show value.(y_2)
    @show value.(x_3)
    @show value.(w_3)
    @show value.(y_3)
    return
end
