abstract Manipulator{StateType, InputType, OutputType} <: DynamicalSystem{StateType, InputType, OutputType}

@generated function dynamics{T, NS, NI, StateType}(robot::Manipulator{StateType}, time, state::State{NS, T}, input::Input{NI, T})
    return quote
        H, C, B = manipulator_dynamics(robot, state)
        H_inv = inv(H)
        tau::Vec{$(length(velocity_type(robot))), T} = B * input - C
        vdot::Vec{$(length(velocity_type(robot))), T} = H_inv * tau
        v = convert(Velocity, state)
        StateType($([:(v[$(i)]) for i in 1:length(velocity_type(robot))]...), $([:(vdot[$(i)]) for i in 1:length(velocity_type(robot))]...))
    end
end
