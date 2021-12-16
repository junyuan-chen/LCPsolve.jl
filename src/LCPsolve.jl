module LCPsolve

using LinearAlgebra
using Printf
using SparseArrays

import Base: show, push!, getindex

export LCP, solve!

include("solver_state_results.jl")
include("solver.jl")

end
