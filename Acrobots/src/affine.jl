immutable LinearSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs} <: DynamicalSystem{StateType, InputType, OutputType}
    A::Mat{NStates, NStates, T}
    B::Mat{NStates, NInputs, T}
    C::Mat{NOutputs, NStates, T}
    D::Mat{NOutputs, NInputs, T}
end

immutable AffineSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs} <: DynamicalSystem{StateType, InputType, OutputType}
    A::Mat{NStates, NStates, T}
    B::Mat{NStates, NInputs, T}
    C::Mat{NOutputs, NStates, T}
    D::Mat{NOutputs, NInputs, T}
    x0::StateType
    u0::InputType
    xd0::StateType
    y0::OutputType
end
AffineSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}(
        A::Mat{NStates, NStates, T},
        B::Mat{NStates, NInputs, T},
        C::Mat{NOutputs, NStates, T},
        D::Mat{NOutputs, NInputs, T},
        x0::StateType,
        u0::InputType,
        xd0::StateType,
        y0::OutputType) = AffineSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}(
            A, B, C, D, x0, u0, xd0, y0)
state_type(sys::AffineSystem) = typeof(sys.x0)

### Helper methods for affine systems. These allow general affine systems to be interpolated using the Interpolations.jl package
call{T}(::Type{Mat{0, 0, T}}, x::Number) = Mat{0,0,T}()
call{N, T}(::Type{Mat{0, N, T}}) = Mat{0, N, T}(tuple([tuple() for i in 1:N]...))
call{N, T}(::Type{Mat{0, N, T}}, x::Number) = Mat{0, N, T}()
call{M, T}(::Type{Mat{M, 0, T}}, x::Number) = Mat{M, 0, T}()
one{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}(
    ::Type{AffineSystem{T, StateType, InputType, OutputType, NStates, NInputs, NOutputs}}) =
        AffineSystem(
    Mat{NStates, NStates, T}(1),
    Mat{NStates, NInputs, T}(1),
    Mat{NOutputs, NStates, T}(1),
    Mat{NOutputs, NInputs, T}(1),
    StateType(1),
    InputType(1),
    StateType(1),
    OutputType(1))
*{T}(x::Real, m::Mat{0, 0, T}) = Mat{0, 0, T}()
*{N, T}(x::Real, m::Mat{0, N, T}) = Mat{0, N, T}()
*{M, T}(x::Real, m::Mat{M, 0, T}) = Mat{M, 0, T}()
*(x::Real, sys::AffineSystem) = AffineSystem(x * sys.A,
    x * sys.B,
    x * sys.C,
    x * sys.D,
    x * sys.x0,
    x * sys.u0,
    x * sys.xd0,
    x * sys.y0)
+{T}(m1::Mat{0, 0, T}, m2::Mat{0, 0, T}) = Mat{0, 0, T}()
+{N, T}(m1::Mat{0, N, T}, m2::Mat{0, N, T}) = Mat{0, N, T}()
+{M, T}(m1::Mat{M, 0, T}, m2::Mat{M, 0, T}) = Mat{M, 0, T}()
+(sys1::AffineSystem, sys2::AffineSystem) = AffineSystem(sys1.A + sys2.A,
    sys1.B + sys2.B,
    sys1.C + sys2.C,
    sys1.D + sys2.D,
    sys1.x0 + sys2.x0,
    sys1.u0 + sys2.u0,
    sys1.xd0 + sys2.xd0,
    sys1.y0 + sys2.y0)



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
        return :(Output(sys.C * state + sys.D * input))
    elseif length(state) > 0
        return :(Output(sys.C * state))
    elseif length(input) > 0
        return :(Output(sys.D * input))
    else
        return :(Output())
    end
end

function output{SysType}(sys::Interpolations.GriddedInterpolation{SysType}, t, state, input)
    output(sys[t], t, state, input)
end

function dynamics{SysType}(sys::Interpolations.GriddedInterpolation{SysType}, t, state, input)
    dynamics(sys[t], t, state, input)
end

output{T, State, Input, Output}(sys::AffineSystem{T, State, Input, Output}, t, state::State, input::Input) = output(convert(LinearSystem, sys), t, state - sys.x0, input - sys.u0) + sys.y0

@generated function dynamics{T, State, Input}(sys::DynamicalSystem{State, Input}, x::AbstractVector{T})
    segment_lengths = [1; length(State); length(Input)]
    breaks = cumsum(segment_lengths)
    :(dynamics(sys, x[1],
        State(x[$(breaks[1]+1):$(breaks[2])]),
        Input(x[$(breaks[2]+1):$(breaks[3])])))
end

@generated function output{T, State, Input}(sys::DynamicalSystem{State, Input}, x::AbstractVector{T})
    segment_lengths = [1; length(State); length(Input)]
    breaks = cumsum(segment_lengths)
    :(output(sys, x[1],
        State(x[$(breaks[1]+1):$(breaks[2])]),
        Input(x[$(breaks[2]+1):$(breaks[3])])))
end

function linearize{StateType, InputType, OutputType}(robot::DynamicalSystem{StateType, InputType, OutputType}, time, state::StateType, input::InputType)
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
call{T <: LQRState}(::Type{T}, x::Number) = LQRState{Float64}()

function lqr{T, State, Input, Output}(sys::AffineSystem{T, State, Input, Output}, Q, R, operating_point)
    K = lqr(Matrix{T}(sys.A), Matrix{T}(sys.B), Q, R)
    u0 = Output(operating_point)
    AffineSystem{T, LQRState, Output, Input, 0, length(Output), length(Input)}(
        Mat{0,0,T}(),
        Mat{0, length(Output), T}(),
        Mat{length(Input), 0, T}(),
        Mat{length(Input), length(Output), T}(-K),
        LQRState{T}(),
        u0,
        LQRState{T}(),
        Input(0))
end

immutable TimeVaryingRiccati
    linearizations
    Qs
    Rs
end

function dynamics(sys::TimeVaryingRiccati, t, state, input)
    Q = sys.Qs[t]
    R = sys.Rs[t]
    affine_sys = sys.linearizations[t]
    A = affine_sys.A
    B = affine_sys.B
    S = state
    -(Q - S * B * inv(R) * B' * S + S * A + A' * S)
end

output(sys::TimeVaryingRiccati, t, state, input) = state

function tvlqr(linearizations, xf, Qs, Rs, Qf, Rf)
    # R = Mat(R)
    # Q = Mat(Q)
    Qf = Mat(Qf)
    Rf = Mat(Rf)
    knots = linearizations.knots[1]
    tf = knots[end]
    sysf = linearizations[tf]
    S = Mat(care(Matrix{Float64}(sysf.A), Matrix{Float64}(sysf.B), Matrix{Float64}(Qf), Matrix{Float64}(Rf)))
    tvriccati = TimeVaryingRiccati(linearizations, Qs, Rs)
    K = inv(Rf) * sysf.B' * S'

    T = Float64
    Output = AcrobotOutput
    Input = AcrobotInput
    controllers = [AffineSystem(
        Mat{0,0,T}(),
        Mat{0, length(Output), T}(),
        Mat{length(Input), 0, T}(),
        Mat{length(Input), length(Output), T}(-K),
        LQRState{T}(),
        Output(xf),
        LQRState{T}(),
        Output(0))]

    for i = length(knots):-1:2
        t = knots[i]
        dt = knots[i-1] - knots[i]
        Sdot = dynamics(tvriccati, t, S, 0)
        S += dt * Sdot

        # TODO: this is duplicating the lookup inside dynamics()
        aff_sys = linearizations[t]
        R = Rs[t]
        K = inv(R) * aff_sys.B' * S'


        push!(controllers, AffineSystem(
            Mat{0,0,T}(),
            Mat{0, length(Output), T}(),
            Mat{length(Input), 0, T}(),
            Mat{length(Input), length(Output), T}(-K),
            LQRState{T}(),
            xf,
            LQRState{T}(),
            aff_sys.u0))
    end
    interpolate((knots,), reverse(controllers), Gridded(Linear()))
end
