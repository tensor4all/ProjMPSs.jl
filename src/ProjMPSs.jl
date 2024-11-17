module ProjMPSs

import OrderedCollections: OrderedSet, OrderedDict
using EllipsisNotation
using LinearAlgebra: norm

import ITensors: ITensors, Index, ITensor, dim
import ITensorMPS: ITensorMPS, AbstractMPS, MPS, MPO, siteinds

import FastMPOContractions as FMPOC
#using Quantics

#import TensorCrossInterpolation as TCI
#import TCIAlgorithms as TCIA

include("util.jl")
include("projector.jl")
include("projmps.jl")
include("blockedmps.jl")
include("patching.jl")
include("contract.jl")
#include("mul.jl")
include("adaptivemul.jl")

end
