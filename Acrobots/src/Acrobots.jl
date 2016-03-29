module Acrobots

using FixedSizeArrays
import FixedSizeArrays: isapprox
using Quaternions
using PyLCM
using PyCall
import Base: convert, *
using ForwardDiff
using Flatten
import ControlSystems: lqr

@pyimport drake as lcmdrake

include("system_types.jl")
include("manipulator.jl")
include("drake_vis.jl")
include("acrobot.jl")
include("affine.jl")

end
