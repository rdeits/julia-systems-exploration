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
