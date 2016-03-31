@make_type DoubleIntegratorState State x xdot
@make_type DoubleIntegratorPosition Position x
@make_type DoubleIntegratorVelocity Velocity x
@make_type DoubleIntegratorInput Input u
@make_type DoubleIntegratorOutput Output x xdot

convert{T, PType <: Position}(::Type{PType}, state::DoubleIntegratorState{T}) = DoubleIntegratorPosition{T}(state.x)
convert{T, VType <: Velocity}(::Type{VType}, state::DoubleIntegratorState{T}) = DoubleIntegratorVelocity{T}(state.xdot)

type DoubleIntegrator{T} <: DynamicalSystem{DoubleIntegratorState, DoubleIntegratorInput, DoubleIntegratorOutput}
    mass::T
end

function dynamics(robot::DoubleIntegrator, t, state, input)
    DoubleIntegratorState(state.xdot, input.u / robot.mass)
end

# state_type{T <: DoubleIntegrator}(bot_type::Type{T}) = DoubleIntegratorState
# velocity_type{T <: DoubleIntegrator}(bot_type::Type{T}) = DoubleIntegratorVelocity
#
#
# function manipulator_dynamics{ParamType, T}(robot::DoubleIntegrator{ParamType}, state::DoubleIntegratorState{T})
#     H = Mat{1, 1, T}(robot.mass)
#     Cv = Vec{2, T}(0)
#     B = Mat{1, 1, T}(1)
# end

double_integrator() = DoubleIntegrator{Float64}(1.0)

output(robot::DoubleIntegrator, t, state, input) = DoubleIntegratorOutput(state)

function viewer_load_msg(robot::DoubleIntegrator)
    msg = lcmdrake.lcmt_viewer_load_robot()
    msg[:num_links] = 1

    geom = lcmdrake.lcmt_viewer_geometry_data()
    geom[:type] = geom[:BOX]
    geom[:position] = [0; 0; 0.5]
    quat = qrotation([0.; 1; 0], 0)
    geom[:quaternion] = [quat.s; quat.v1; quat.v2; quat.v3]
    geom[:color] = [0.2; 0.2; 0.8; 0.6]
    geom[:string_data] = ""
    geom[:float_data] = [1; 1; 1]
    geom[:num_float_data] = 3
    link_msg = lcmdrake.lcmt_viewer_link_data()
    link_msg[:name] = "double_integrator"
    link_msg[:robot_num] = 1
    link_msg[:num_geom] = 1
    push!(link_msg["geom"], geom)

    push!(msg["link"], link_msg)
    msg
end

function viewer_draw_msg(robot::DoubleIntegrator, position::DoubleIntegratorPosition)
    msg = lcmdrake.lcmt_viewer_draw()
    msg[:num_links] = 1
    msg[:link_name] = ["double_integrator"]
    msg[:robot_num] = [1]
    msg[:position] = Any[[position.x; 0; 0]]
    msg[:quaternion] = Any[[1; 0; 0; 0]]
    msg
end
