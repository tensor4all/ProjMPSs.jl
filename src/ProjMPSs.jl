module ProjMPSs

import OrderedCollections: OrderedSet, OrderedDict
using EllipsisNotation
using LinearAlgebra: norm

import ITensors: ITensors, MPS, MPO, siteinds, Index, AbstractMPS, ITensor, dim

import FastMPOContractions as FMPOC
using Quantics

#import TensorCrossInterpolation as TCI
#import TCIAlgorithms as TCIA

include("util.jl")
include("projector.jl")
include("projmps.jl")
include("blockedmps.jl")
#include("conversion.jl")
include("patching.jl")
include("contract.jl")
include("mul.jl")
include("adaptivemul.jl")

end
