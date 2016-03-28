type LinearSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}
    A::Mat{NStates, NStates, T}
    B::Mat{NStates, NInputs, T}
    C::Mat{NOutputs, NStates, T}
    D::Mat{NOutputs, NInputs, T}
end

(*){M, N, T}(A::Mat{M, N, T}, x::FixedVectorNoTuple{N, T}) = [sum([A[i,j] * x[j] for j in 1:N]) for i in 1:M]

function dynamics{T, StateType, InputType}(sys::LinearSystem{T, StateType, InputType}, t, state::StateType, input::InputType)
    StateType{T}(sys.A * state + sys.B * input)
end

function output{T, State, Input, Output}(sys::LinearSystem{T, State, Input, Output}, t, state::State, input::Input)
    Output{T}(sys.C * state + sys.D * input)
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

@make_type LQRState State x

function lqr{T, State, Input, Output}(sys::LinearSystem{T, State, Input, Output}, Q, R)
    K = lqr(Matrix{T}(sys.A), Matrix{T}(sys.B), Q, R)
    LinearSystem{T, LQRState, Output, Input, 1, length(Output), length(Input)}(
        Mat{1,1,T}(0),
        Mat{1, length(Output), T}(0),
        Mat{length(Input), 1, T}(0),
        Mat{length(Input), length(Output), T}(-K))
end
