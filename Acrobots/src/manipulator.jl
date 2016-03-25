abstract Manipulator{StateType, InputType, OutputType}

function dynamics{State, Input}(robot::Manipulator, time, state::State, input::Input)
    H, C, B = manipulator_dynamics(robot, state)
    H_inv = inv(H)
    tau = B * destructure([input]) - C
    vdot = H_inv * tau
    State(convert(Velocity, state)..., vdot...)
end
