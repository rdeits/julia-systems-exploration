# workspace()
module Acrobots

using FixedSizeArrays
using Quaternions
using PyLCM
using PyCall
import Base: convert
using ForwardDiff
using Flatten

@pyimport drake as lcmdrake

macro make_type(typename, parent, fields...)
    ex = :(immutable $(esc(typename)){T} <: $(parent){$(length(fields)), T} end)
    push!(ex.args, Expr(:block))
    for field in fields
        push!(ex.args[3].args, :($(field)::T))
    end
    constructor_expr = quote function $(esc(typename))(a::NTuple{$(length(fields)), T})
            new{T}()
            end
        end
    for i = 1:length(fields)
        push!(constructor_expr.args[2].args[2].args[2].args, :(a[$(i)]))
    end
    push!(ex.args[3].args, constructor_expr)
    return ex
end

type State{PositionType, VelocityType}
    position::PositionType
    velocity::VelocityType
end
Base.length{PositionType, VelocityType}(T::State{PositionType, VelocityType}) = length(PositionType) + length(VelocityType)
abstract Position{N, T} <: FixedVectorNoTuple{N, T}
abstract Velocity{N, T} <: FixedVectorNoTuple{N, T}
abstract Input{N, T} <: FixedVectorNoTuple{N, T}
abstract Output{N, T} <: FixedVectorNoTuple{N, T}

convert{PosType <: Position, VelType <: Velocity}(::Type{PosType}, state::State{PosType, VelType}) = state.position
convert{PosType <: Position, VelType <: Velocity}(::Type{VelType}, state::State{PosType, VelType}) = state.velocity

abstract Manipulator{StateType, InputType, OutputType}

function dynamics{State, Input}(robot::Manipulator, time, state::State, input::Input)
    H, C, B = manipulator_dynamics(robot, state)
    H_inv = inv(H)
    tau = B * destructure([input]) - C
    vdot = H_inv * tau
    State(state.velocity, vec(vdot))
end

type DrakeVisualizer
    robot::Manipulator
    lcm::PyLCM.PyLCMWrapper

    function DrakeVisualizer(robot::Manipulator, lcm::PyLCM.PyLCMWrapper)
        vis = new(robot, lcm)
        load_robot(vis)
        vis
    end
end

DrakeVisualizer(robot::Manipulator) = DrakeVisualizer(robot, LCM())

function load_robot(vis::DrakeVisualizer)
    load_msg = viewer_load_msg(vis.robot)
    publish(vis.lcm, "DRAKE_VIEWER_LOAD_ROBOT", load_msg)
end

function draw{T}(vis::DrakeVisualizer, state::State{T})
    publish(vis.lcm, "DRAKE_VIEWER_DRAW", viewer_draw_msg(vis.robot, state))
end

@make_type AcrobotPosition Position theta1 theta2
@make_type AcrobotVelocity Velocity theta1 theta2
@make_type AcrobotInput Input tau
@make_type AcrobotOutput Output theta1 theta2 theta1dot theta2dot
typealias AcrobotState{T} State{AcrobotPosition{T}, AcrobotVelocity{T}}

type AcrobotLink{T}
    length::T
    mass::T
    damping::T
    length_to_CoM::T
    inertia::T
end

type Acrobot{T} <: Manipulator{AcrobotState, AcrobotInput, AcrobotOutput}
    links::NTuple{2, AcrobotLink{T}}
    gravity::T
end

function manipulator_dynamics{ParamType, T}(robot::Acrobot{ParamType}, state::AcrobotState{T})
    inertias_about_joint = [robot.links[i].inertia + robot.links[i].mass * robot.links[i].length_to_CoM^2 for i in 1:2]
    m2l1lc2 = robot.links[2].mass * robot.links[1].length * robot.links[2].length_to_CoM

    c = cos(state.position)
    s = sin(state.position)
    s12 = sin(state.position.theta1 + state.position.theta2)

    h12 = inertias_about_joint[2] + m2l1lc2 * c[2]
    H = [inertias_about_joint[1] + inertias_about_joint[2] + robot.links[2].mass * robot.links[1].length^2 + 2 * m2l1lc2 * c[2]   h12;
        h12                                                                                                                       inertias_about_joint[2]]
    C = [-2 * m2l1lc2 * s[2] * state.velocity.theta2     -m2l1lc2 * s[2] * state.velocity.theta2;
         m2l1lc2 * s[2] * state.velocity.theta1          0]
    G = robot.gravity * [robot.links[1].mass * robot.links[1].length_to_CoM * s[1] + robot.links[2].mass * (robot.links[1].length * s[1] + robot.links[2].length_to_CoM * s12);
                         robot.links[2].mass * robot.links[2].length_to_CoM * s12]

    Cv = Float64[dot(AcrobotVelocity{T}(vec(C[i,:])), state.velocity) + G[i] + robot.links[i].damping * state.velocity[i] for i = 1:2]

    B = Float64[0; 1]

    H, Cv, B
end

acrobot() = Acrobot((AcrobotLink(1.0, 1.0, 0.1, 0.5, 0.083), AcrobotLink(2.0, 1.0, 0.1, 1.0, 0.33)), 9.81)

function output{T}(robot::Acrobot, time, state::AcrobotState{T}, input::AcrobotInput{T})
    AcrobotOutput{T}(state.position[1], state.position[2], state.velocity[1], state.velocity[2])
end


function viewer_data(link::Acrobots.AcrobotLink, name)
    msg = lcmdrake.lcmt_viewer_link_data()
    msg[:name] = name
    msg[:robot_num] = 1
    msg[:num_geom] = 1

    geom = lcmdrake.lcmt_viewer_geometry_data()
    geom[:type] = geom[:BOX]
    geom[:position] = [0; 0; -link.length / 2]
    quat = qrotation([0.; 1; 0], pi/2)
    geom[:quaternion] = [quat.s; quat.v1; quat.v2; quat.v3]
    geom[:color] = [0.2; 0.2; 0.8; 0.6]
    geom[:string_data] = ""
    geom[:float_data] = [link.length; 0.1; 0.1]
    geom[:num_float_data] = 3
    push!(msg["geom"], geom)

    msg
end

function viewer_load_msg(robot::Acrobot)
    msg = lcmdrake.lcmt_viewer_load_robot()
    msg[:num_links] = 2
    for i in 1:2
        push!(msg["link"], viewer_data(robot.links[i], "link$(i)"))
    end
    msg
end


rotmat(theta) = [cos(theta) -sin(theta); sin(theta) cos(theta)]

function viewer_draw_msg{T}(robot::Acrobot, state::AcrobotState{T})
    msg = lcmdrake.lcmt_viewer_draw()
    msg[:num_links] = 2
    msg[:link_name] = ["link1"; "link2"]
    msg[:robot_num] = [1; 1]

    p0 = [0; 0]
    p1 = rotmat(state.position.theta1) * [0; -robot.links[1].length]

    msg[:position] = Any[[0; 0; 0], [p1[1]; 0; p1[2]]]
    quats = [qrotation([0; -1; 0], state.position.theta1);
            qrotation([0.; -1; 0], state.position.theta1 + state.position.theta2)]
    msg[:quaternion] = Any[[q.s; q.v1; q.v2; q.v3] for q in quats]
    msg
end

type LinearSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}
    A::Mat{NStates, NStates, T}
    B::Mat{NStates, NInputs, T}
    C::Mat{NOutputs, NStates, T}
    D::Mat{NOutputs, NInputs, T}
end

function dynamics{T, StateType, InputType}(sys::LinearSystem{T, StateType, InputType}, t, state::StateType, input::InputType)
    StateType{T}(sys.A * [state.position; state.velocity] + sys.B * input)
end

function linearize{StateType, InputType, OutputType}(robot::Manipulator{StateType, InputType, OutputType}, time, state::StateType, input::InputType)
    segment_lengths = [1; length(state.position) + length(state.velocity); length(input)]
    breaks = cumsum(segment_lengths)
    @show StateType
    @show InputType
    function wrapped_dynamics(x)
        to_vector(dynamics(robot, x[1],
            x[(breaks[1]+1):breaks[2]],
            x[(breaks[2]+1):breaks[3]]))
    end
    function wrapped_output(x)
        to_vector(output(robot, x[1],
            x[(breaks[1]+1):breaks[2]],
            x[(breaks[2]+1):breaks[3]]))
    end

    x = vcat(map(x -> convert(Vector{Float64}, x), ([time], state.position, state.velocity, input))...)
    AB = ForwardDiff.jacobian(wrapped_dynamics, x)
    CD = FowardDiff.jacobian(wrapped_output, x)
    LinearSystem{Float64, StateType, InputType, OutputType, length(StateType), length(InputType), length(OutputType)}(
        AB[:,1:length(StateType)], AB[:,(length(StateType+1):end)],
        CD[:,1:length(StateType)], CD[:,(length(StateType+1):end)])
end

end
