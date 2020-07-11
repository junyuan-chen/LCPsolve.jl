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
