type LinearSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}
    A::Mat{NStates, NStates, T}
    B::Mat{NStates, NInputs, T}
    C::Mat{NOutputs, NStates, T}
    D::Mat{NOutputs, NInputs, T}
end

@generated function dynamics{T, StateType, InputType}(sys::LinearSystem{T, StateType, InputType}, t, state::StateType, input::InputType)
    if length(state) > 0 && length(input) > 0
        return :(StateType{T}(sys.A * state + sys.B * input))
    elseif length(state) > 0
        return :(StateType{T}(sys.A * state))
    elseif length(input) > 0
        return :(StateType{T}(sys.B * input))
    else
        return :(StateType{T}())
    end
end

@generated function output{T, State, Input, Output}(sys::LinearSystem{T, State, Input, Output}, t, state::State, input::Input)
    if length(state) > 0 && length(input) > 0
        return :(Output{T}(sys.C * state + sys.D * input))
    elseif length(state) > 0
        return :(Output{T}(sys.C * state))
    elseif length(input) > 0
        return :(Output{T}(sys.D * input))
    else
        return :(Output{T}())
    end
end

function linearize{StateType, InputType, OutputType}(robot::Manipulator{StateType, InputType, OutputType}, time, state::StateType, input::InputType)
    segment_lengths = [1; length(state); length(input)]
    breaks = cumsum(segment_lengths)
    function wrapped_dynamics(x)
        to_vector(dynamics(robot, x[1],
            StateType(x[(breaks[1]+1):breaks[2]]),
            InputType(x[(breaks[2]+1):breaks[3]])))
    end
    function wrapped_output(x)
        to_vector(output(robot, x[1],
            StateType(x[(breaks[1]+1):breaks[2]]),
            InputType(x[(breaks[2]+1):breaks[3]])))
    end

    x = vcat(map(x -> convert(Vector{Float64}, x), ([time], state, input))...)
    AB = ForwardDiff.jacobian(wrapped_dynamics, x)
    CD = ForwardDiff.jacobian(wrapped_output, x)
    LinearSystem{Float64, StateType, InputType, OutputType, length(StateType), length(InputType), length(OutputType)}(
        AB[:,1+(1:length(StateType))], AB[:,((length(StateType)+2):end)],
        CD[:,1+(1:length(StateType))], CD[:,((length(StateType)+2):end)])
end

@make_type LQRState State

function lqr{T, State, Input, Output}(sys::LinearSystem{T, State, Input, Output}, Q, R)
    K = lqr(Matrix{T}(sys.A), Matrix{T}(sys.B), Q, R)
    LinearSystem{T, LQRState, Output, Input, 0, length(Output), length(Input)}(
        Mat{0,0,T}(),
        Mat{0, length(Output), T}(tuple([tuple() for i in 1:length(Output)]...)),
        Mat{length(Input), 0, T}(),
        Mat{length(Input), length(Output), T}(-K))
end
