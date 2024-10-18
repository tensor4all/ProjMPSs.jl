module ProjMPSs

import OrderedCollections: OrderedSet
using EllipsisNotation

using ITensors

import FastMPOContractions as FMPOC

using Quantics
#import Quantics: MatrixMultiplier, ElementwiseMultiplier, preprocess, postprocess, asMPS

include("util.jl")
include("projector.jl")
include("blockstructure.jl")
include("projmps.jl")

include("blockedmps.jl")

include("contract.jl")
include("mul.jl")

end
