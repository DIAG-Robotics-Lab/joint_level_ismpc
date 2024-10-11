#pragma once

#include "Controller.hpp"

#include <dart/dart.hpp>
//#include <dart/utils/utils.hpp>
#include <dart/gui/osg/osg.hpp>

#include <dart/gui/gui.hpp>


class WorldNode : public dart::gui::osg::WorldNode
{
public:
  // Constructor
  WorldNode(dart::simulation::WorldPtr world, dart::dynamics::SkeletonPtr robot, std::shared_ptr<Controller> controller);

  // Documentation inherited
  void customPreStep() override;

protected:
  std::shared_ptr<Controller> mController;
  Eigen::Vector3d mExternalForce;
  int mForceDuration;
  dart::simulation::WorldPtr mWorld;
  dart::dynamics::SkeletonPtr mRobot;

  dart::dynamics::SimpleFramePtr mTarget; //TODO understand if needed
};

