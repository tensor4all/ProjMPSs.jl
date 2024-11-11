module ProjMPSs

import OrderedCollections: OrderedSet, OrderedDict
using EllipsisNotation

import ITensors: ITensors, MPS, MPO, siteinds, Index, AbstractMPS, ITensor
using LinearAlgebra: norm
import FastMPOContractions as FMPOC

using Quantics

include("util.jl")
include("projector.jl")
include("projmps.jl")

include("blockedmps.jl")

include("patching.jl")

include("contract.jl")
include("mul.jl")
include("adaptivemul.jl")

end
