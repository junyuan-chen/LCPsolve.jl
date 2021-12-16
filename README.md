# LCPsolve.jl

*A solver for linear complementarity problems*

[![CI-stable][CI-stable-img]][CI-stable-url]
[![codecov][codecov-img]][codecov-url]
[![PkgEval][pkgeval-img]][pkgeval-url]

[CI-stable-img]: https://github.com/junyuan-chen/LCPsolve.jl/workflows/CI-stable/badge.svg
[CI-stable-url]: https://github.com/junyuan-chen/LCPsolve.jl/actions?query=workflow%3ACI-stable

[codecov-img]: https://codecov.io/gh/junyuan-chen/LCPsolve.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/junyuan-chen/LCPsolve.jl

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/L/LCPsolve.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/L/LCPsolve.html

[LCPsolve.jl](https://github.com/junyuan-chen/LCPsolve.jl)
provides a Julia implementation of the
[Matlab routine](https://www.mathworks.com/matlabcentral/fileexchange/20952-lcp-mcp-solver-newton-based)
written by Yuval Tassa.
The solver is particularly useful when the problem to be solved is ill-conditioned.
This is often the case,
for example,
when the linear system arises from a discretization of
a Hamilton-Jacobi-Bellman variational inequality
in the process of value function iteration.
Illustrative applications in economics for solving optimal stopping problems can be found
[here](https://benjaminmoll.com/codes/).

Please see [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl)
for solvers suitable for nonlinear problems.

## Usage

An object of type `LCP` is used to specify the problem.
Passing the object to `solve!` yields the results,
which is returned in an object of type `SolverResults`.
The solution is stored in the field `sol`.

```julia
using LCPsolve, SparseArrays
n = 10000
M = spdiagm(0=>[9;fill(17,n-2);9], -1=>fill(-8,n-1), 1=>fill(-8,n-1))
q = -log.(collect(LinRange(0.05,5,n)))
result = solve!(LCP(M,q))
```

#### Setting up the Problem

The complementarity problem can be specified
by constructing an object of type `LCP`.
By default, calling `LCP(M, q)` creates a problem for finding a vector `x` such that
for a square matrix `M` and a vector `q`
the following three equations hold simultaneously:

```math
            x >= 0;
       Mx + q >= 0;
    x'(Mx + q) = 0.
```

It is advisable to construct `M` as a sparse matrix of type `SparseMatrixCSC`.

By providing `LCP` the keyword arguments `l` and `u`,
which represent constraints as vectors,
one can solve the more general problem
in which each element of `x` indexed by `i`
satisfies the following relations:

```math
    l[i] < x[i] < u[i]   =>   (Mx)[i] + q[i] = 0;
    (Mx)[i] + q[i] < 0   =>   x[i] = u[i];
    (Mx)[i] + q[i] > 0   =>   x[i] = l[i].
```

#### Solving the Problem

To solve the problem, pass the `LCP` object to `solve!`.
The method `solve!` accepts an optional argument `x0`
as the initial starting point for the solver.
Keyword arguments can be passed for adjusting the behavior of the solver.
For details, please use the help mode in REPL.

An object of type `SolverResults` will be returned by `solve!`.
The solution is stored in the field `sol`.
The field `converged` indicates whether convergence has been reached.

#### Applications in Optimal Stopping Problems

The `examples` folder contains more illustrations for the usage of the solver.

## References

Fischer, A. (1995). A Newton-type method for positive-semidefinite linear complementarity problems. Journal of Optimization Theory and Applications, 86(3), 585-608.

Bazaraa, M. S., Sherali, H. D., & Shetty, C. M. (2013). Nonlinear programming: Theory and algorithms. John Wiley & Sons.

Tassa, Y. (2008). LCP / MCP solver (Newton-based). MATLAB Central File Exchange. Retrieved from https://www.mathworks.com/matlabcentral/fileexchange/20952-lcp-mcp-solver-newton-based.
