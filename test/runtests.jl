import ProjMPSs as PM
using Random
using Test

include("_util.jl")

include("projector_tests.jl")
include("blockstructure_tests.jl")
include("projmps_tests.jl")
include("blockedmps_tests.jl")
include("contract_tests.jl")
include("mul_tests.jl")
