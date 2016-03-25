@make_type AcrobotState State theta1 theta2 theta1dot theta2dot
@make_type AcrobotPosition Position theta1 theta2
@make_type AcrobotVelocity Velocity theta1 theta2
@make_type AcrobotInput Input tau
@make_type AcrobotOutput Output theta1 theta2 theta1dot theta2dot
# typealias AcrobotState{T} State{AcrobotPosition{T}, AcrobotVelocity{T}}

convert{T, PType <: Position}(::Type{PType}, state::AcrobotState{T}) = AcrobotPosition{T}(state.theta1, state.theta2)
convert{T, VType <: Velocity}(::Type{VType}, state::AcrobotState{T}) = AcrobotVelocity{T}(state.theta1dot, state.theta2dot)

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

    position = convert(AcrobotPosition, state)
    velocity = convert(AcrobotVelocity, state)

    c = cos(position)
    s = sin(position)
    s12 = sin(position.theta1 + position.theta2)

    h12 = inertias_about_joint[2] + m2l1lc2 * c[2]
    H = [inertias_about_joint[1] + inertias_about_joint[2] + robot.links[2].mass * robot.links[1].length^2 + 2 * m2l1lc2 * c[2]   h12;
        h12                                                                                                                       inertias_about_joint[2]]
    C = [-2 * m2l1lc2 * s[2] * velocity.theta2     -m2l1lc2 * s[2] * velocity.theta2;
         m2l1lc2 * s[2] * velocity.theta1          0]
    G = robot.gravity * [robot.links[1].mass * robot.links[1].length_to_CoM * s[1] + robot.links[2].mass * (robot.links[1].length * s[1] + robot.links[2].length_to_CoM * s12);
                         robot.links[2].mass * robot.links[2].length_to_CoM * s12]

    Cv = T[dot(AcrobotVelocity{T}(vec(C[i,:])), velocity) + G[i] + robot.links[i].damping * velocity[i] for i = 1:2]

    B = Float64[0; 1]

    H, Cv, B
end

acrobot() = Acrobot((AcrobotLink(1.0, 1.0, 0.1, 0.5, 0.083), AcrobotLink(2.0, 1.0, 0.1, 1.0, 0.33)), 9.81)

function output{T}(robot::Acrobot, time, state::AcrobotState{T}, input::AcrobotInput{T})
    AcrobotOutput{T}(state...)
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
    p1 = rotmat(state.theta1) * [0; -robot.links[1].length]

    msg[:position] = Any[[0; 0; 0], [p1[1]; 0; p1[2]]]
    quats = [qrotation([0; -1; 0], state.theta1);
            qrotation([0.; -1; 0], state.theta1 + state.theta2)]
    msg[:quaternion] = Any[[q.s; q.v1; q.v2; q.v3] for q in quats]
    msg
end
