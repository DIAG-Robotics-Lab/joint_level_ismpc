#include "WorldNode.hpp"

#include "dart/external/imgui/imgui.h"

WorldNode::WorldNode(const dart::simulation::WorldPtr world,
                     const dart::dynamics::SkeletonPtr robot,
                     const std::shared_ptr<Controller> controller)
    : mWorld(world), mRobot(robot), mController(controller), dart::gui::osg::WorldNode(world),
      mExternalForce(Eigen::Vector3d::Zero()), mForceDuration(0.0) {

    //mController.reset(new Controller(hrp4, world));
    //mController->setInitialConfiguration();

}

void WorldNode::customPreStep() {
    mController->update();
}

