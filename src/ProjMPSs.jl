module ProjMPSs

import OrderedCollections: OrderedSet, OrderedDict
using EllipsisNotation
using LinearAlgebra: LinearAlgebra

import ITensors: ITensors, Index, ITensor, dim, inds, qr, commoninds, uniqueinds
import ITensorMPS: ITensorMPS, AbstractMPS, MPS, MPO, siteinds, findsites
import ITensors.TagSets: hastag, hastags

import FastMPOContractions as FMPOC

include("util.jl")
include("projector.jl")
include("projmps.jl")
include("blockedmps.jl")
include("patching.jl")
include("contract.jl")
include("adaptivemul.jl")

# Only for backward compatibility
include("automul.jl")

end
