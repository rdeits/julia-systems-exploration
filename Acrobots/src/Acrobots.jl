module Acrobots

using FixedSizeArrays
using DrakeVisualizer
import DrakeVisualizer: draw
using AffineTransforms
using Quaternions
using ForwardDiff
using Interpolations
import GeometryTypes: HyperRectangle, Vec
import FixedSizeArrays: isapprox
import ControlSystems: lqr, care
import Base: convert, one, *, +, call

include("system_types.jl")
include("manipulator.jl")
include("drake_vis.jl")
include("acrobot.jl")
include("double_integrator.jl")
include("affine.jl")
include("SLQ.jl")

end
