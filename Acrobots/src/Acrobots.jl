module Acrobots

using FixedSizeArrays
import FixedSizeArrays: isapprox
using Quaternions
using PyLCM
using PyCall
import Base: convert, one, *, +, call
using ForwardDiff
import ControlSystems: lqr, care
using Interpolations

@pyimport drake as lcmdrake

include("system_types.jl")
include("manipulator.jl")
include("drake_vis.jl")
include("acrobot.jl")
include("affine.jl")

end
