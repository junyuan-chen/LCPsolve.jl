module LCPsolver

using LinearAlgebra
using SparseArrays
using Printf

import Base.show,
       Base.push!,
       Base.getindex,
       Base.setindex!

export LCP, solve!

include("solver_state_results.jl")
include("solver.jl")

end
