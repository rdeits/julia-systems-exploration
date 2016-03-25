type DrakeVisualizer
    robot::Manipulator
    lcm::PyLCM.PyLCMWrapper

    function DrakeVisualizer(robot::Manipulator, lcm::PyLCM.PyLCMWrapper)
        vis = new(robot, lcm)
        load_robot(vis)
        vis
    end
end

DrakeVisualizer(robot::Manipulator) = DrakeVisualizer(robot, LCM())

function load_robot(vis::DrakeVisualizer)
    load_msg = viewer_load_msg(vis.robot)
    publish(vis.lcm, "DRAKE_VIEWER_LOAD_ROBOT", load_msg)
end

function draw{T}(vis::DrakeVisualizer, state::State{T})
    publish(vis.lcm, "DRAKE_VIEWER_DRAW", viewer_draw_msg(vis.robot, state))
end
