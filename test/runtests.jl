using Test
using LCPsolve
using LinearAlgebra
using SparseArrays

tests = ["solve.jl",
         "examples.jl"]

println("Running tests:")

for test in tests
    include(test)
end
