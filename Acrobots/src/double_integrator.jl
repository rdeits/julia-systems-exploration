@make_type DoubleIntegratorState State x xdot
@make_type DoubleIntegratorPosition Position x
@make_type DoubleIntegratorVelocity Velocity x
@make_type DoubleIntegratorInput Input u
@make_type DoubleIntegratorOutput Output x xdot

convert{T, PType <: Position}(::Type{PType}, state::DoubleIntegratorState{T}) = DoubleIntegratorPosition{T}(state.x)
convert{T, VType <: Velocity}(::Type{VType}, state::DoubleIntegratorState{T}) = DoubleIntegratorVelocity{T}(state.xdot)

type DoubleIntegrator{T} <: DynamicalSystem{DoubleIntegratorState, DoubleIntegratorInput, DoubleIntegratorOutput}
    mass::T
end

function dynamics(robot::DoubleIntegrator, t, state, input)
    DoubleIntegratorState(state.xdot, input.u / robot.mass)
end

double_integrator() = DoubleIntegrator{Float64}(1.0)

output(robot::DoubleIntegrator, t, state, input) = DoubleIntegratorOutput(state)

function visualizer_load(robot::DoubleIntegrator)
    geometry = HyperRectangle(Vec(-0.5,-0.5,0), Vec(1.,1,1))
    model = DrakeVisualizer.load(geometry)
    return model
end

link_origins(robot::DoubleIntegrator, position::DoubleIntegratorPosition) = [tformtranslate([position.x; 0; 0])]
