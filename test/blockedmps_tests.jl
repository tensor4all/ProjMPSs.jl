using Test

using ITensors

using Quantics: Quantics
import ProjMPSs: Projector, project, ProjMPS, BlockedMPS

@testset "blockedmps.jl" begin
    @testset "two blocks" begin
        N = 3
        sitesx = [Index(2, "x=$n") for n in 1:N]
        sitesy = [Index(2, "y=$n") for n in 1:N]

        sites = collect(collect.(zip(sitesx, sitesy)))

        Ψ = MPS(collect(_random_mpo(sites)))

        prjΨ = ProjMPS(Ψ)

        prjΨ1 = project(prjΨ, Dict(sitesx[1] => 1))
        prjΨ2 = project(prjΨ, Dict(sitesx[1] => 2))

        @test_throws ErrorException BlockedMPS([prjΨ, prjΨ1])
        @test_throws ErrorException BlockedMPS([prjΨ1, prjΨ1])

        Ψreconst = BlockedMPS(prjΨ1) + BlockedMPS(prjΨ2)
        @test length(Ψreconst.data) == 2
        @test MPS(Ψreconst) ≈ Ψ
    end
end
