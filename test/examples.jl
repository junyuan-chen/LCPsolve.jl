# Test the scripts in the examples folder

@testset "car" begin
    include("../examples/car.jl")
    # Compare results with the Matlab code
    @test isapprox(H.mV[1,:], [-201.8163489871447, -192.4728306882469])
    @test isapprox(H.mC[end,:], [0.238122006231259, 0.247896044803184])
end

@testset "Hopenhayn" begin
    include("../examples/Hopenhayn.jl")
    # Compare results with the Matlab code
    @test isapprox(dot(f.vN, f.vG.*g.dz), 0.768574232318884)
    @test findfirst(x->x!=0.0, f.vV) == 564
end