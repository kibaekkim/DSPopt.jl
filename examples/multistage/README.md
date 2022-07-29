# Multistage Stochastic Optimization

The example is described in page 270 of the following book:

- Birge and Louveaux. "Introduction to Stochastic Programming, Second Edition."

## Problem Description

Consider a 3-stage scenario tree with a single root node at stage 1 and two nodes in stage 2, each of which has two children nodes in stage 3 (thus four leaf nodes). Let $\mathcal{N}$ denote the set of nodes, $\mathcal{L} \subset \mathcal{N}$ the set of leaf nodes in the scenario tree, and $\alpha(n)$ the ancestor of node $n$ (with $n = 1$ denoting the root node). The air conditioning production planning problem from Birge & Louveaux can be modeled as a 3-stage stochastic integer program:

$$\eqalign{
\min & h_{1}^{\top} {\left\lbrack \matrix{x_{1} \cr w_{1} \cr y_{1}} \right\rbrack} + \sum_{n \in \mathcal{N} \setminus \\{\\{1\\}, \{\mathcal{L}\}\\}} 0.5 h_{2}^{\top} {\left\lbrack \matrix{x_{n} \cr w_{n} \cr y_{n}} \right\rbrack} + \sum_{n \in \{\mathcal{L}\}} 0.25 h_{3}^{\top} {\left\lbrack \matrix{x_{n} \cr w_{n} \cr y_{n}} \right\rbrack} \\
\text{s.t.} \ & {\left\lbrack \matrix{1 & 1 & -1} \right\rbrack} {\left\lbrack \matrix{x_{1} \cr w_{1} \cr y_{1}} \right\rbrack} = d_{1}, \\
& {\left\lbrack \matrix{1 & 1 & 1 & -1} \right\rbrack} {\left\lbrack \matrix{ y_{\alpha(n)} \cr x_{n} \cr w_{n} \cr y_{n}} \right\rbrack} = d_{n} \quad \forall n \in \mathcal{N} \setminus \\{1\\}, \\
& x_{n}, w_{n}, y_{n} \in H_n,
}$$ 

where
$h_{1}^{\top} = h_{2}^{\top} = {\left\lbrack \matrix{1 & 3 & 0.5} \right\rbrack}, \ h_{3}^{\top} = {\left\lbrack \matrix{1 & 3 & 0} \right\rbrack}, \ d^{\top} = {\left\lbrack \matrix{1 & 1 & 3 & 1 & 3 & 1 & 3} \right\rbrack}, \ \text{and} \ H_{n} = \\{(x_{n},w_{n},y_{n}) \in \mathbb{Z}^{3}_{+} \ | \ x_n \leq 2 \\}.$

We optimize this model using a scenario-tree decomposition via the following 4 methods:

### A deterministic integer programming model

Let $x_t[k]$ be the regular-time production at time $t$ in scenario $k$, $w_t[k]$ be the overtime production at time $t$ in scenario $k$, and $y_t[k]$ be the number of units stored at time $t$ in scenario $k$.

Solution for JuMP model

- Optimal objective value: 6.25
- $x_1^* = 2, x_2^* = [1, 2]^{\top}, x_3^* = [0, 2, 1, 2]^{\top}$
- $w_1^* = 0, w_2^* = [0, 0]^{\top}, w_3^* = [0, 0, 0, 1]^{\top}$
- $y_1^* = 1, y_2^* = [1, 0]^{\top}, y_3^* = [0, 0, 0, 0]^{\top}$

### A three-stage integer programming model

Let $x[n,t]$ be the regular-time production at node $n$ at time period $t$, $w[n,t]$ be the overtime production at node $n$ at time period $t$, and $y[n,t]$ be the number of units stored at node $n$ at time period $t$.

3. Solution for DSP model __with__ non-anticipativity.
- Optimal objective value: 7.25
- $x^* = [2, 2, 2, 1, 1, 2, 2]^{\top}$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^{\top}$
- $y^* = [1, 2, 0, 2, 0, 1, 0]^{\top}$

4. Solution for DSP model __without__ non-anticipativity. 
- Optimal objective value: 6.25
- $x^* = [2, 1, 2, 0, 2, 1, 2]^{\top}$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^{\top}$
- $y^* = [1, 1, 0, 0, 0, 0, 0]^{\top}$

For the methods __without__ non-anticipativity, we assume that both regular-time labor and overtime labor can be decided after the uncertainty is realized. For the methods __with__ non-anticipativity, we assume that regular-time labor must be decided before the uncertainty is realized, and that overtime labor is a corrective measure that can be taken once the uncertainty is realized - therefore, the variable encoding regular-time labor, $x$, is included in the non-anticipativity constraints, while the variable encoding overtime labor, $w$, is not. 

The JuMP model can be found in "jump_test.jl", and the DSP model can be found in "three_stages.jl". 
