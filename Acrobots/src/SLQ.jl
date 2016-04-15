function simulate_and_linearize{State, Input, Output}(robot::DynamicalSystem{State, Input, Output},
    initial_state,
    controller,
    ts::Range,
    ticks_per_draw=1,
    vis=Visualizer(robot))

    dt = step(ts)
    x = initial_state
    u = Input(0.0)
    linearizations = AffineSystem{Float64, State, Input, Output, length(State), length(Input), length(Output)}[]

    @assert length(state_type(controller)) == 0
    controller_state = state_type(controller)(0.0)
    for (i, t) in enumerate(ts)
        y = output(robot, t, x, u)
        u = output(controller, t, controller_state, y)
        xdot = dynamics(robot, t, x, u)
        x += xdot * dt
        push!(linearizations, linearize(robot, t, x, u))

        if mod(i, ticks_per_draw) == 0
            draw(vis, x)
        end
    end

    return linearizations
end

function quadratize_costs{State, Input, Output}(robot::DynamicalSystem{State, Input, Output}, ts::Range, linearizations, state_cost, input_cost)
    Q_generator = ForwardDiff.hessian(x -> state_cost(State(x)), ForwardDiff.AllResults)
    R_generator = ForwardDiff.hessian(x -> input_cost(Input(x)), ForwardDiff.AllResults)

    Qs = Mat{length(State), length(State), Float64}[]
    Rs = Mat{length(Input), length(Input), Float64}[]
    qs = Vec{length(State), Float64}[]
    rs = Vec{length(Input), Float64}[]

    for (i, t) in enumerate(ts)
        state = linearizations[i].x0
        input = linearizations[i].u0
        Q, Q_results = Q_generator(convert(Vector, state))
        push!(Qs, Q)
        push!(qs, ForwardDiff.gradient(Q_results))
        R, R_results = R_generator(convert(Vector, input))
        push!(Rs, R)
        push!(rs, ForwardDiff.gradient(R_results))
    end

    Qs, qs, Rs, rs
end

function search_direction(robot, ts::Range, linearizations, Qs, qs, Rs, rs, S_final, xdes)
    dt = step(ts)
    Ps = Array{Mat{2, 2, Float64}}(length(ts))
    ps = Array{Vec{2, Float64}}(length(ts))
    Ks = Array{Mat{1, 2, Float64}}(length(ts))
    ls = Array{Vec{1, Float64}}(length(ts))
    As = [Mat(eye(2)) + dt * sys.A for sys in linearizations]
    Bs = [dt * sys.B for sys in linearizations]

    Ps[end] = S_final
    ps[end] = S_final * (linearizations[end].x0 - xdes)

    for j = (length(ts)-1):-1:1
        g = rs[j] + Bs[j]' * ps[j+1]
        G = Bs[j]' * Ps[j+1] * As[j]
        H = Rs[j] + Bs[j]' * Ps[j+1] * Bs[j]
        Hi = inv(H)
        Ks[j] = -Hi * G
        ls[j] = -Hi * g

        Ps[j] = Qs[j] + As[j]' * Ps[j+1] * As[j] + Ks[j]' * H * Ks[j] + Ks[j]' * G + G' * Ks[j]
        ps[j] = qs[j] + As[j]' * ps[j+1] + Ks[j]' * H * ls[j] + Ks[j]' * g + G' * ls[j]
    end
    Ks[end] = 0
    ls[end] = 0

    Ks, ls
end

function update_controller{State, Input, Output}(robot::DynamicalSystem{State, Input, Output}, ts::Range, linearizations, Ks, ls, alpha=0.5)
    controllers = AffineSystem{Float64, LQRState{Float64}, Output{Float64},
        Input{Float64}, 0, length(Output), length(Input)}[AffineSystem(
        Mat{0,0,Float64}(),
        Mat{0,2,Float64}(),
        Mat{1,0,Float64}(),
        Mat{1,2,Float64}(Ks[j]),
        LQRState{Float64}(),
        Output{Float64}(linearizations[j].x0),
        LQRState{Float64}(),
        Input{Float64}(linearizations[j].u0) + Input{Float64}(alpha * ls[j]))
        for j in 1:length(ts)]

    interpolate((ts,), controllers, Gridded(Constant()))
end
