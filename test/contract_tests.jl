using Test
import ProjMPSs: Projector, ProjMPS, projcontract

@testset "contract.jl" begin
    @testset "contract (xk-y-z)" begin
        R = 3
        sitesx = [Index(2, "Qubit,x=$n") for n in 1:R]
        sitesk = [Index(2, "Qubit,k=$n") for n in 1:R]
        sitesy = [Index(2, "Qubit,y=$n") for n in 1:R]
        sitesz = [Index(2, "Qubit,z=$n") for n in 1:R]

        sitesa = collect(collect.(zip(sitesx, sitesk, sitesy)))
        sitesb = collect(collect.(zip(sitesy, sitesz)))

        p1 = ProjMPS(_random_mpo(sitesa), Projector(Dict(sitesx[1] => 1)))
        p2 = ProjMPS(_random_mpo(sitesb), Projector(Dict(sitesz[1] => 1)))

        p12 = ITensors.contract(p1, p2; alg="naive")

        @test p12.projector == Projector(Dict(sitesx[1] => 1, sitesz[1] => 1))

        proj_subset = Projector(Dict(sitesx[1] => 1, sitesz[1] => 1, sitesk[1] => 1))
        p12_2 = projcontract(p1, p2, proj_subset; alg="naive")

        @test p12_2.projector == proj_subset
    end
end
