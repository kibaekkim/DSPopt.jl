# Multistage Stochastic Optimization

The example is described in page 270 of the following book:

- Birge and Louveaux. "Introduction to Stochastic Programming, Second Edition."

## Problem Description

__Scenario Tree__: Consider a 3-stage scenario tree with a single root node at stage 1 and two nodes in stage 2, each of which has two children nodes in stage 3 (thus four leaf nodes). Let $\mathcal{N}$ denote the set of nodes, $\mathcal{L} \subset \mathcal{N}$ the set of leaf nodes in the scenario tree, and $\alpha(n)$ the ancestor of node $n$ (with $n = 1$ denoting the root node).
<br/><br/>
__Scenario Lattice__: Consider a 3-stage scenario lattice with a single root node at stage 1 with two children nodes in stage 2, and two leaf nodes in stage 3 that are each connected to both of the nodes in stage 2. Let $\mathcal{N}$ denote the set of nodes, $\mathcal{L} \subset \mathcal{N}$ the set of leaf nodes in the scenario lattice (with $n = 1$ denoting the root node). 
<br/><br/>
The air conditioning production planning problem from Birge & Louveaux can be modeled as a stochastic integer program using both a scenario tree or a scenario lattice, either by nodal or scenario decomposition:

#### Nodal Decomposition using a Scenario Tree: 
<center>Let $x_{n}^{t}$ be the regular-time production at node $n$ at time period $t$, $w_{n}^{t}$ be the overtime production at node $n$ at time period $t$, and $y_{n}^{t}$ be the number of units stored at node $n$ at time period $t$. We assume there are no units stored before $t=0$, so $y_{1}^{-1} = 0$.<center> 

$$\eqalign{
\min & h_{1}^{\top} {\left\lbrack \matrix{x_{1}^{0} \cr w_{1}^{0} \cr y_{1}^{0}} \right\rbrack} + \sum_{n \in \mathcal{N} \setminus \\{\\{1\\}, \{\mathcal{L}\}\\}} 0.5 h_{2}^{\top} {\left\lbrack \matrix{x_{n}^{1} \cr w_{n}^{1} \cr y_{n}^{1}} \right\rbrack} + \sum_{n \in \{\mathcal{L}\}} 0.25 h_{3}^{\top} {\left\lbrack \matrix{x_{n}^{2} \cr w_{n}^{2}} \right\rbrack} \\
\text{s.t.} \ & {\left\lbrack \matrix{1 & 1 & 1 & -1} \right\rbrack} {\left\lbrack \matrix{ y_{n}^{t-1} \cr x_{n}^{t} \cr w_{n}^{t} \cr y_{n}^{t}} \right\rbrack} = d_{n} \quad \forall n \in \mathcal{N}, \\
& x_{n}^{t}, w_{n}^{t}, y_{n}^{t} \in H_n \quad \forall n \in \mathcal{N}, 
}$$ 

<center>where $h_{1}^{\top} = h_{2}^{\top} = {\left\lbrack \matrix{1 & 3 & 0.5} \right\rbrack}, \ h_{3}^{\top} = {\left\lbrack \matrix{1 & 3} \right\rbrack}, \ d^{\top} = {\left\lbrack \matrix{1 & 1 & 3 & 1 & 3 & 1 & 3} \right\rbrack}, \ \text{and} \ H_{n} = \\{(x_{n},w_{n},y_{n}) \in \mathbb{Z}^{3}_{+} \ | \ x_n \leq 2 \\}.$<center>

<center>We can incorporate non-anticipativity into the model endogenously by having $t$ correspond to the stage of node $n$ for all decision variables. We can also model non-anticipativity explicitly by having:
- $t$ denote the stage of node $n$ for variables $x$ and $w$.
- Decision variables $y_{n}^{t}$ AND $y_{n}^{t-1}$ for each node $n$, where $t$ denotes the stage of node $n$. 
- Including the non-anticipativity constraints: $y_{\lfloor \frac{n}{2} \rfloor}^{t-1} = y_{n}^{t-1}$.<center>

#### Nodal Decomposition using a Scenario Lattice: 
  
<center>Let $x_{n}^{t}$ be the regular-time production at node $n$ at time period $t$, $w_{n}^{t}$ be the overtime production at node $n$ at time period $t$, and $y_{n}^{t}$ be the number of units stored at node $n$ at time period $t$. We assume there are no units stored before $t=0$.<center> 

$$\eqalign{
\min & h_{1}^{\top} {\left\lbrack \matrix{x_{1}^{0} \cr w_{1}^{0} \cr y_{1}^{0}} \right\rbrack} + \sum_{n \in \mathcal{N} \setminus \\{\\{1\\}, \{\mathcal{L}\}\\}} 0.5 h_{2}^{\top} {\left\lbrack \matrix{x_{n}^{1} \cr w_{n}^{1} \cr y_{n}^{1}} \right\rbrack} + \sum_{n \in \{\mathcal{L}\}} 0.5 h_{3}^{\top} {\left\lbrack \matrix{x_{n}^{2} \cr w_{n}^{2}} \right\rbrack} \\
\text{s.t.} \ & {\left\lbrack \matrix{1 & 1 & 1 & -1} \right\rbrack} {\left\lbrack \matrix{ y_{n}^{t-i} \cr x_{n}^{t} \cr w_{n}^{t} \cr y_{n}^{t}} \right\rbrack} = d_{n} \quad \forall n \in \mathcal{N}, i \in \\{1,...,t\\}, \\
& x_{n}^{t}, w_{n}^{t}, y_{n}^{t} \in H_n \quad \forall n \in \mathcal{N}. 
}$$ 

<center>Let $\mathcal{N^{even}}$ be the set of even-numbered nodes (similarly with $\mathcal{N^{odd}}$). We can, again, incorporate non-anticipativity endogenously into the model, or include the constraints: $y_{n}^{t} = y_{n+i}^{t-1} \ \forall n \in \\{\mathcal{N^{even}} \setminus \mathcal{L}\\}, i = \\{2,3\\}$, and $y_{n}^{t} = y_{n+i}^{t} \ \forall n \in \\{\mathcal{N^{odd}} \setminus \mathcal{L}\\}, i = \\{1,2\\}$.<center>

#### Scenario Decomposition using a Scenario Tree:

<center>Let $x_{t}^{k}$ be the regular-time production at time $t$ in scenario $k$, $w_{t}^{k}$ be the overtime production at time $t$ in scenario $k$, and $y_{t}^{k}$ be the number of units stored at time $t$ in scenario $k$. Since we assume there are no units stored before $t=0$, $y_{\alpha(1)} = 0$.<center> 

$$\eqalign{
\min & h_{1}^{\top} {\left\lbrack \matrix{x_{1}^{1} \cr w_{1}^{1} \cr y_{1}^{1}} \right\rbrack} + \sum_{n \in \mathcal{N} \setminus \\{\\{1\\}, \{\mathcal{L}\}\\}, \ k = 1,2} 0.5 h_{2}^{\top} {\left\lbrack \matrix{x_{n}^{k} \cr w_{n}^{k} \cr y_{n}^{k}} \right\rbrack} + \sum_{n \in \{\mathcal{L}\}, \ k = 1, 2, 3, 4} 0.25 h_{3}^{\top} {\left\lbrack \matrix{x_{n}^{k} \cr w_{n}^{k}} \right\rbrack} \\
\text{s.t.} \ & {\left\lbrack \matrix{1 & 1 & 1 & -1} \right\rbrack} {\left\lbrack \matrix{ y_{\alpha(n)}^{k} \cr x_{n}^{k} \cr w_{n}^{k} \cr y_{n}^{k}} \right\rbrack} = d_{n} \quad \forall n \in \mathcal{N}, \\
& x_{n}^{k}, w_{n}^{k}, y_{n}^{k} \in H_n \quad \forall n \in \mathcal{N}.
}$$ 

<center>Similarly with the nodal decomposition, we can incorporate non-anticipativity into the model endogenously or explicitly.<center> 

## File Descriptions

The files in this folder include: 

#### A deterministic integer programming model using a scenario tree (scenario decomposition)

Solution for JuMP model

- Optimal objective value: 6.25
- $x_1^* = 2, x_2^* = [1, 2]^{\top}, x_3^* = [0, 2, 1, 2]^{\top}$
- $w_1^* = 0, w_2^* = [0, 0]^{\top}, w_3^* = [0, 0, 0, 1]^{\top}$
- $y_1^* = 1, y_2^* = [1, 0]^{\top}, y_3^* = [0, 0, 0, 0]^{\top}$

#### A deterministic integer programming model using a scenario tree (nodal decomposition)

Solution for JuMP model

- Optimal objective value: 6.25
- $x^* = [2, 1, 2, 0, 2, 1, 2]^{\top}$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^{\top}$
- $y^* = [1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0]^{\top}$

#### A three-stage integer programming model using a scenario tree (nodal decomposition)

Solution for DSP model 

- Optimal objective value: 6.25
- $x^* = [2, 1, 2, 0, 2, 1, 2]^{\top}$
- $w^* = [0, 0, 0, 0, 0, 0, 1]^{\top}$
- $y^* = [1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0]^{\top}$

#### A three-stage integer programming model using a scenario lattice (nodal decomposition)

Solution for DSP model 

- Optimal objective value: 6.5
- $x^* = [2, 0, 2, 1, 2]^{\top}$
- $w^* = [0, 0, 0, 0, 1]^{\top}$
- $y^* = [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]^{\top}$

Solution for JuMP model

- Optimal objective value: 6.5
- $x^* = [2, 0, 2, 1, 2]^{\top}$
- $w^* = [0, 0, 0, 0, 1]^{\top}$
- $y^* = [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]^{\top}$

The JuMP models can be found in [``` scen_decomp.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/scen_decomp.jl) (scen decomp w/ scen tree), [``` node_decomp.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/node_decomp.jl) (nodal decomp w/ scen tree), and [``` jump_lattice.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/jump_lattice.jl) (nodal decomp w/ scen lattice). The DSP models can be found in [``` DSP_tree.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/DSP_tree.jl) (nodal decomp w/ scen tree) and [``` DSP_lattice.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/DSP_lattice.jl) (nodal decomp w/ scen lattice). 

## Notes

- Each model (except [``` scen_decomp.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/scen_decomp.jl)) can be extended to more than 3 stages by changing the ``` T = ``` parameter in each file to the number of stages of your choosing. 
- In [``` DSP_tree.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/DSP_tree.jl), the ``` num_scenarios ``` parameter corresponds to the number of scenarios represented in the scenario tree. However, in [``` DSP_lattice.jl ```](https://github.com/kibaekkim/DSPopt.jl/blob/ra/multistage/examples/multistage/DSP_lattice.jl), the same parameter corresponds to the number of nodes in the scenario lattice. 
- For the DSP models, each variable must have an upper bound. 
