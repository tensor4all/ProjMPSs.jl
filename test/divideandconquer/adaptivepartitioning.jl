@testitem "PatchOrdering" begin
    po = TCIA.PatchOrdering([4, 3, 2, 1])
    @test TCIA.maskactiveindices(po, 2) == [1, 1, 0, 0]
    @test TCIA.maskactiveindices(po, 1) == [1, 1, 1, 0]
    @test TCIA.fullindices(po, [[1]], [[2], [3], [4]]) == [[2], [3], [4], [1]]
    @test TCIA.fullindices(po, [[1], [2]], [[3], [4]]) == [[3], [4], [2], [1]]
end

@testitem "2D fermi gk" for _flipordering in [false, true]
    using Distributed
    using Test

    # Define the maximum number of worker processes.
    const MAX_WORKERS = 4

    # Add worker processes if necessary.
    addprocs(max(0, MAX_WORKERS - nworkers()))

    @everywhere using TensorCrossInterpolation
    @everywhere import TensorCrossInterpolation as TCI
    @everywhere import TCIAlgorithms as TCIA
    @everywhere using ITensors
    @everywhere ITensors.disable_warn_order()


    Random.seed!(1234)

    ek(kx, ky) = 2 * cos(kx) + 2 * cos(ky) - 1.0

    function gk(kx, ky, β)
        iv = im * π / β
        return 1 / (iv - ek(kx, ky))
    end

    R = 20
    grid = DiscretizedGrid{2}(R, (0.0, 0.0), (2π, 2π))
    localdims = fill(4, R)

    β = 20.0
    flipper = _flipordering ? x -> reverse(x) : x -> x
    f = x -> gk(originalcoordinate(grid, QuanticsInd{2}.(flipper(x)))..., β)

    tol = 1e-5

    pordering = TCIA.PatchOrdering(flipper(collect(1:R)))

    creator = TCIA.TCI2PatchCreator(
        ComplexF64,
        f,
        localdims;
        maxbonddim=100,
        rtol=tol,
        verbosity=1,
        ntry=10,
    )

    tree = TCIA.adaptivepatches(creator, pordering; verbosity=1, maxnleaves=1000)
    @show tree

    _evaluate(x, idx) = FMPOC.evaluate(x, [[i] for i in idx])

    for _ = 1:100
        pivot = rand(1:4, R)
        error_func = x -> abs(f(x) - _evaluate(tree, x))
        pivot = TCI.optfirstpivot(error_func, localdims, pivot)
        @test isapprox(_evaluate(tree, pivot), f(pivot); atol=10 * creator.atol)
    end

end

nothing
