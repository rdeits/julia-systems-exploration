type DrakeVisualizer
    robot
    lcm::PyLCM.PyLCMWrapper

    function DrakeVisualizer(robot, lcm::PyLCM.PyLCMWrapper)
        vis = new(robot, lcm)
        load_robot(vis)
        vis
    end
end

DrakeVisualizer(robot) = DrakeVisualizer(robot, LCM())

function load_robot(vis::DrakeVisualizer)
    load_msg = viewer_load_msg(vis.robot)
    publish(vis.lcm, "DRAKE_VIEWER_LOAD_ROBOT", load_msg)
end

function draw{T}(vis::DrakeVisualizer, position::Position{T})
    publish(vis.lcm, "DRAKE_VIEWER_DRAW", viewer_draw_msg(vis.robot, position))
end

function draw{T}(vis::DrakeVisualizer, state::State{T})
    draw(vis, convert(Position, state))
end
