module Acrobots

using FixedSizeArrays
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
include("linear.jl")

end
