using Test
import ProjMPSs: Projector, project, ProjMPS, projcontract, _add_directsum
import FastMPOContractions as FMPOC
import Quantics: asMPO

@testset "contract.jl" begin
    @testset "contract (xk-y-z)" begin
        R = 3
        sitesx = [Index(2, "Qubit,x=$n") for n in 1:R]
        sitesk = [Index(2, "Qubit,k=$n") for n in 1:R]
        sitesy = [Index(2, "Qubit,y=$n") for n in 1:R]
        sitesz = [Index(2, "Qubit,z=$n") for n in 1:R]

        sitesa = collect(collect.(zip(sitesx, sitesk, sitesy)))
        sitesb = collect(collect.(zip(sitesy, sitesz)))

        p1 = project(ProjMPS(_random_mpo(sitesa)), Projector(Dict(sitesx[1] => 1)))
        p2 = project(ProjMPS(_random_mpo(sitesb)), Projector(Dict(sitesz[1] => 1)))

        p12 = ITensors.contract(p1, p2; alg="naive")

        @test p12.projector == Projector(Dict(sitesx[1] => 1, sitesz[1] => 1))

        proj_subset = Projector(Dict(sitesx[1] => 1, sitesz[1] => 1, sitesk[1] => 1))
        p12_2 = projcontract(p1, p2, proj_subset; alg="naive")

        @test p12_2.projector == proj_subset
    end

    @testset "contract (2x2)" begin
        R = 10
        cutoff = 1e-25
        linkdims = 2

        sitesx = [Index(2, "Qubit,x=$n") for n in 1:R]
        sitesy = [Index(2, "Qubit,y=$n") for n in 1:R]
        sitesz = [Index(2, "Qubit,z=$n") for n in 1:R]

        sitesa = collect(collect.(zip(sitesx, sitesy)))
        sitesb = collect(collect.(zip(sitesy, sitesz)))

        a = _random_mpo(sitesa; linkdims=linkdims)
        b = _random_mpo(sitesb; linkdims=linkdims)

        proj_a = [
            project(ProjMPS(a), Projector(Dict(sitesx[1] => i, sitesy[1] => j))) for
            i in 1:2, j in 1:2
        ]
        proj_b = [
            project(ProjMPS(b), Projector(Dict(sitesy[1] => i, sitesz[1] => j))) for
            i in 1:2, j in 1:2
        ]

        for x in [1, 2], y in [1, 2]
            res = projcontract(
                vec(proj_a),
                vec(proj_b),
                Projector(Dict(sitesx[1] => x, sitesz[1] => y));
                cutoff=cutoff,
            )
            ref = reduce(
                _add_directsum,
                [
                    FMPOC.contract_mpo_mpo(
                        asMPO(proj_a[x, k].data), asMPO(proj_b[k, y].data); alg="naive"
                    ) for k in 1:2
                ],
            )
            @test res.data ≈ MPS(collect(ref.data))
        end

        ab = contract(BlockedMPS(vec(proj_a)), BlockedMPS(vec(proj_b)); alg="naive")
        @test MPS(ab) ≈ MPS(collect(contract(a, b; alg="naive")))
    end
end
