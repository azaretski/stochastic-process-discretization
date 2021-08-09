# stochastic-process-discretization
Gauss-Hermite quadrature and discretization of independent AR(1)-processes to a joint discrete-state Markov chain
- `gauss_hermite` computes n Gauss-Hermite quadrature nodes with the corresponding weights;
- `get_MC_ind` computes indices of a discretized Markov chain based on a sequence of uniform draws;
- `joint_MC` discretizes k independent AR(1) processes to discrete-state Markov chains and constructs a joint Markov chain and transition matrix.

`Julia` versions of `joint_MC` and `get_MC_ind` are in `discretized_MC.jl`.
