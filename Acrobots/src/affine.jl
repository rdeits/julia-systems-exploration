immutable LinearSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}
    A::Mat{NStates, NStates, T}
    B::Mat{NStates, NInputs, T}
    C::Mat{NOutputs, NStates, T}
    D::Mat{NOutputs, NInputs, T}
end

immutable AffineSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}
    A::Mat{NStates, NStates, T}
    B::Mat{NStates, NInputs, T}
    C::Mat{NOutputs, NStates, T}
    D::Mat{NOutputs, NInputs, T}
    x0::StateType
    u0::InputType
    xd0::StateType
    y0::OutputType
end

convert{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}(::Type{LinearSystem}, sys::AffineSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}) = LinearSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}(sys.A, sys.B, sys.C, sys.D)

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

function dynamics{T, StateType, InputType}(sys::AffineSystem{T, StateType, InputType}, t, state::StateType, input::InputType)
    return dynamics(convert(LinearSystem, sys), t, state - sys.x0, input - sys.u0) + sys.xd0
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

output{T, State, Input, Output}(sys::AffineSystem{T, State, Input, Output}, t, state::State, input::Input) = output(convert(LinearSystem, sys), t, state - sys.x0, input - sys.u0) + sys.y0

@generated function dynamics{T, State, Input}(sys::Manipulator{State, Input}, x::AbstractVector{T})
    segment_lengths = [1; length(State); length(Input)]
    breaks = cumsum(segment_lengths)
    :(dynamics(sys, x[1],
        State(x[$(breaks[1]+1):$(breaks[2])]),
        Input(x[$(breaks[2]+1):$(breaks[3])])))
end

@generated function output{T, State, Input}(sys::Manipulator{State, Input}, x::AbstractVector{T})
    segment_lengths = [1; length(State); length(Input)]
    breaks = cumsum(segment_lengths)
    :(output(sys, x[1],
        State(x[$(breaks[1]+1):$(breaks[2])]),
        Input(x[$(breaks[2]+1):$(breaks[3])])))
end

function linearize{StateType, InputType, OutputType}(robot::Manipulator{StateType, InputType, OutputType}, time, state::StateType, input::InputType)
    x = vcat(map(x -> convert(Vector{Float64}, x), ([time], state, input))...)
    AB, dynamics_results = ForwardDiff.jacobian(x -> dynamics(robot, x), x, ForwardDiff.AllResults)
    CD, output_results = ForwardDiff.jacobian(x -> output(robot, x), x, ForwardDiff.AllResults)
    AffineSystem{Float64, StateType, InputType, OutputType, length(StateType), length(InputType), length(OutputType)}(
        AB[:,1+(1:length(StateType))], AB[:,((length(StateType)+2):end)],
        CD[:,1+(1:length(StateType))], CD[:,((length(StateType)+2):end)],
        state,
        input,
        StateType(ForwardDiff.value(dynamics_results)),
        OutputType(ForwardDiff.value(output_results)))
end

@make_type LQRState State

function lqr{T, State, Input, Output}(sys::AffineSystem{T, State, Input, Output}, Q, R)
    K = lqr(Matrix{T}(sys.A), Matrix{T}(sys.B), Q, R)
    AffineSystem{T, LQRState, Output, Input, 0, length(Output), length(Input)}(
        Mat{0,0,T}(),
        Mat{0, length(Output), T}(tuple([tuple() for i in 1:length(Output)]...)),
        Mat{length(Input), 0, T}(),
        Mat{length(Input), length(Output), T}(-K),
        LQRState{T}(),
        sys.C * sys.x0 + sys.D * sys.u0,
        LQRState{T}(),
        Output(0))
end
