using Test
import ProjMPSs: ProjMPSs, Projector, project, ProjMPS, adaptive_patching, BlockedMPS
import FastMPOContractions as FMPOC
using ITensors
using Random

@testset "patching.jl" begin
    @testset "adaptive_patching" begin
        Random.seed!(1234)

        R = 3
        sitesx = [Index(2, "Qubit,x=$n") for n in 1:R]
        sitesy = [Index(2, "Qubit,y=$n") for n in 1:R]

        sites = collect(collect.(zip(sitesx, sitesy)))

        prjmps = ProjMPS(_random_mpo(sites; linkdims=20))

        sites_ = collect(Iterators.flatten(sites))
        bmps = BlockedMPS(adaptive_patching(prjmps, sites_; maxdim=10, cutoff=1e-25))

        @test length(values((bmps))) > 1

        @test MPS(bmps) ≈ MPS(prjmps) rtol = 1e-15
        #MPS(bmps) ≈ MPS(prjmps)
    end
end
