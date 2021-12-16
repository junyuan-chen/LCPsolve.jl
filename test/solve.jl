# Randomly generate a strictly diagonally dominant matrix M and a vector q
function randn_Mq(n::Int, sparse::Bool)
    M = sparse ? sprandn(n, n, 5/n) : randn(n, n)
    M[diagind(M)] .= sum(abs.(M), dims=2)[:].+ones(n)
    q = randn(n)
    return M, q
end

@testset "diagonally dominant" begin
    n = 1000
    i = 1
    for sp in [true, false],
        lu in [(zeros(n), fill(Inf,n)), (fill(-Inf,n), zeros(n)), (fill(-Inf,n), fill(Inf,n))]
        # The case with both l and u being finite is not guaranteed to converge
        M, q = randn_Mq(n, true)
        println("Test case ", i, ":")
        i += 1
        lcp = LCP(M, q, l=lu[1], u=lu[2])
        s = solve!(lcp, show_trace=true, max_iter=50)
        @test s.converged
    end
end

@testset "state results" begin
    M = [2.3032505278991975 -0.7864190378835538 -0.09394002534918648;
        1.4271438213161987 2.9612966358546906 0.06405614545993135;
        -2.2064646539344093 -1.3261919640113726 6.687356756945661]
    q = [0.08329207979985552, 0.4644471809076496, -0.45215372730490144]
    lcp = LCP(M, q, l=zeros(3), u=[10.0,Inf,Inf])
    s = solve!(lcp, store_trace=true, max_iter=50)
    # Compre results with Matlab
    # LCP(M, q, zeros(3,1), reshape([Inf; Inf; Inf], 3,1), zeros(3,1), true)
    @test s.sol ≈ [0.0, 0.0, 0.06761320320569002] atol=1e-6
    @test sprint(show, s.trace[1]) == "iter =  1, ψ = 1e-01, r = 0.9, μ = 2e-04\n"
    @test sprint(show, s.trace) == """
        iter =  1, ψ = 1e-01, r = 0.9, μ = 2e-04
        iter =  2, ψ = 6e-04, r = 1.0, μ = 4e-05
        iter =  3, ψ = 9e-06, r = 1.0, μ = 8e-06
        iter =  4, ψ = 8e-09, r = 1.0, μ = 0e+00
        iter =  5, ψ = 7e-15, r = 1.0, μ = 0e+00"""
    @test sprint(show, s) == """
        Results of Solving LCP
         * Convergence: true
         * iter =  5, ψ = 7e-15, r = 1.0, μ = 0e+00"""
end
