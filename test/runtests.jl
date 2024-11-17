import ProjMPSs as PM
using Random
using ITensors
using ITensorMPS
using Test

include("_util.jl")

include("projector_tests.jl")
include("projmps_tests.jl")
include("blockedmps_tests.jl")
include("contract_tests.jl")
include("patching_tests.jl")
