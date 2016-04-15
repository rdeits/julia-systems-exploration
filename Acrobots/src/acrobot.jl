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

velocity_type{T <: Acrobot}(bot_type::Type{T}) = AcrobotVelocity

function manipulator_dynamics{ParamType, T}(robot::Acrobot{ParamType}, state::AcrobotState{T})
    inertias_about_joint = [robot.links[i].inertia + robot.links[i].mass * robot.links[i].length_to_CoM^2 for i in 1:2]
    m2l1lc2 = robot.links[2].mass * robot.links[1].length * robot.links[2].length_to_CoM

    position = convert(AcrobotPosition{T}, state)
    velocity = convert(AcrobotVelocity{T}, state)

    c = cos(position)
    s = sin(position)
    s12 = sin(position.theta1 + position.theta2)

    h12 = inertias_about_joint[2] + m2l1lc2 * c[2]
    H = Mat{2, 2, T}((inertias_about_joint[1] + inertias_about_joint[2] + robot.links[2].mass * robot.links[1].length^2 + 2 * m2l1lc2 * c[2], h12),
            (h12, inertias_about_joint[2]))

    C = Mat{2, 2, T}((-2 * m2l1lc2 * s[2] * velocity.theta2, m2l1lc2 * s[2] * velocity.theta1),
                     (-m2l1lc2 * s[2] * velocity.theta2, 0))

    G = robot.gravity * Vec{2, T}(robot.links[1].mass * robot.links[1].length_to_CoM * s[1] + robot.links[2].mass * (robot.links[1].length * s[1] + robot.links[2].length_to_CoM * s12),
                         robot.links[2].mass * robot.links[2].length_to_CoM * s12)
    damping = Vec{2, T}(robot.links[1].damping, robot.links[2].damping)
    v = Vec{2, T}(velocity)
    Cv::Vec{2, T} = (C * v) + G + damping .* v

    B = Mat{1, 2, T}(0., 1.)'

    H, Cv, B
end

acrobot() = Acrobot((AcrobotLink(1.0, 1.0, 0.1, 0.5, 0.083), AcrobotLink(2.0, 1.0, 0.1, 1.0, 0.33)), 9.81)

function output{T}(robot::Acrobot, time, state::AcrobotState{T}, input::AcrobotInput{T})
    AcrobotOutput{T}(state...)
end

function visualizer_load(robot::Acrobot)
    links = DrakeVisualizer.Link[]
    for i in 1:2
        geometry = HyperRectangle(Vec(0.,0,0), Vec(robot.links[i].length, 0.1, 0.1))
        data = DrakeVisualizer.GeometryData(geometry, tformtranslate([0; -0.05; -0.05]))
        push!(links, DrakeVisualizer.Link([data], "link$(i)"))
    end
    return DrakeVisualizer.load(links)
end

function link_origins(robot::Acrobot, position::AcrobotPosition)
    transforms = Array{AffineTransform{Float64, 3}}(2)
    transforms[1] = tformrotate([0; position.theta1 + pi/2; 0])
    transforms[2] = transforms[1] * tformtranslate([robot.links[1].length; 0; 0]) * tformrotate([0; position.theta2; 0])
    transforms
end
