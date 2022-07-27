Consider a 3-stage scenario tree with a single root node at stage 1 and two nodes in stage 2, each of which has two children nodes in stage 3 (thus four leaf nodes). We solve the air conditioning production planning problem from Birge & Louveaux using a scenario-tree decomposition via the following 4 methods:

### A single-stage integer programming model
Let $x[n]$ be the regular-time production at node $n$, $w[n]$ be the overtime production at node $n$, and $y[n]$ be the number of units stored at node $n$.

1. Solution for JuMP model __with__ non-anticipativity. 
- Optimal objective value: 7.25
- $x^* = [2, 2, 2, 1, 1, 2, 2]^T$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^T$
- $y^* = [1, 2, 0, 2, 0, 1, 0]^T$

2. Solution for JuMP model __without__ non-anticipativity. 
- Optimal objective value: 6.25
- $x^* = [2, 1, 2, 0, 2, 1, 2]^T$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^T$
- $y^* = [1, 1, 0, 0, 0, 0, 0]^T$

### A three-stage integer programming model
Let $x[n,t]$ be the regular-time production at node $n$ at time period $t$, $w[n,t]$ be the overtime production at node $n$ at time period $t$, and $y[n,t]$ be the number of units stored at node $n$ at time period $t$.

2. Solution for DSP model __with__ non-anticipativity.
- Optimal objective value: 7.25
- $x^* = [2, 2, 2, 1, 1, 2, 2]^T$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^T$
- $y^* = [1, 2, 0, 2, 0, 1, 0]^T$

3. Solution for DSP model __without__ non-anticipativity. 
- Optimal objective value: 6.25
- $x^* = [2, 1, 2, 0, 2, 1, 2]^T$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^T$
- $y^* = [1, 1, 0, 0, 0, 0, 0]^T$

We assume that overtime labor is a corrective measure that can be taken once the uncertainty is realized - therefore, the variable encoding overtime labor, $w$, is not included in the non-anticipativity constraints. 

The JuMP model can be found in "jump_test.jl", and the DSP model can be found in "three_stage.jl". 
