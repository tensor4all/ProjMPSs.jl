using Test
import ProjMPSs: BlockStructure, Projector

@testset "BlockStructure" begin
    @testset "constructor" begin
        inds = [Index(2, "n=$n") for n in 1:3]

        let
            p1 = Projector()
            @test_throws ErrorException BlockStructure([p1, p1])
        end

        let
            p1 = Projector(Dict(inds[1] => 1))
            p2 = Projector(Dict(inds[1] => 2))
            bs = BlockStructure([p1, p2])
            @test length(bs) == 2
        end
    end

    #==
    @testset "iterator" begin
        sitedims = [[2], [2], [2]]

        p1 = Projector([[1], [0], [0]], sitedims)
        p2 = Projector([[2], [0], [0]], sitedims)
        bs = BlockStructure([p1, p2])
        @test all(collect(bs) .== [p1, p2])
    end
    ==#
end
