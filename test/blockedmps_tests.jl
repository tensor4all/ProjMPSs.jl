using Test

using ITensors

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
        @test Ψreconst[1] == prjΨ1
        @test Ψreconst[2] == prjΨ2
        @test MPS(Ψreconst) ≈ Ψ
        @test ITensors.norm(Ψreconst) ≈ ITensors.norm(MPS(Ψreconst))
    end

    @testset "two blocks (general key)" begin
        N = 3
        sitesx = [Index(2, "x=$n") for n in 1:N]
        sitesy = [Index(2, "y=$n") for n in 1:N]

        sites = collect(collect.(zip(sitesx, sitesy)))

        Ψ = MPS(collect(_random_mpo(sites)))

        prjΨ = ProjMPS(Ψ)

        prjΨ1 = project(prjΨ, Dict(sitesx[1] => 1))
        prjΨ2 = project(prjΨ, Dict(sitesx[1] => 2))

        a = BlockedMPS(prjΨ1)
        b = BlockedMPS(prjΨ2)

        @test MPS(2 * a) ≈ 2 * MPS(a) rtol = 1e-13
        @test MPS(a * 2) ≈ 2 * MPS(a) rtol = 1e-13
        @test MPS((a + b) + 2 * (b + a)) ≈ 3 * Ψ rtol = 1e-13
        @test MPS((a + b) + 2 * (b + a)) ≈ 3 * Ψ rtol = 1e-13
    end
end
