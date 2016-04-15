type Visualizer{T <: DynamicalSystem}
    robot::T
    model::DrakeVisualizer.VisualizerModel
end

Visualizer(robot) = Visualizer(robot, visualizer_load(robot))

function draw{T}(vis::Visualizer, position::Position{T})
    draw(vis.model, link_origins(vis.robot, position))
end

function draw{T}(vis::Visualizer, state::State{T})
    draw(vis, convert(Position, state))
end
