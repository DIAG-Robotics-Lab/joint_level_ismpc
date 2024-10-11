
#include "WorldNode.hpp"
#include "parameters.hpp"

#include <dart/gui/osg/osg.hpp>
#include <dart/utils/urdf/urdf.hpp>
#include <string>

#include <labrob_qpsolvers/qpsolvers.hpp>

#include <iostream>
#include <algorithm>
#include <filesystem>



int main()
{

  // Create a world
  dart::simulation::WorldPtr world(new dart::simulation::World);

  const std::string urdf_filename = std::filesystem::absolute(models_path+"/hrp4.urdf");
  const std::string ground_filename = std::filesystem::absolute(models_path+"/ground.urdf");
  // const std::string urdf_filename = realpath("../../urdf/hrp4_payload.urdf", nullptr);

  dart::utils::DartLoader::Options options;
  options.mDefaultInertia = dart::dynamics::Inertia(1e-3, Eigen::Vector3d::Zero(), 1e-4*Eigen::Matrix3d::Identity());
  dart::utils::DartLoader urdfLoader(options);

  // Load ground and HRP4 robot and add them to the world
  auto ground = urdfLoader.parseSkeleton(ground_filename);
  auto hrp4 = urdfLoader.parseSkeleton(urdf_filename);
  world->addSkeleton(ground);
  world->addSkeleton(hrp4);

  // PINOCCHIO MODEL



  // Set gravity of the world
  world->setGravity(Eigen::Vector3d(0.0, 0.0, -9.81));
  world->setTimeStep(worldTimeStep);

  // set joint actuator type and force limits
  double forceLimit = 100;

  for (int i = 0; i < hrp4->getNumJoints(); i++) {
    size_t dim = hrp4->getJoint(i)->getNumDofs();
    if(dim==6) {
      hrp4->getJoint(i)->setActuatorType(dart::dynamics::Joint::PASSIVE);
    }
    if(dim==1) {
//      hrp4->getJoint(i)->setActuatorType(dart::dynamics::Joint::VELOCITY);
      hrp4->getJoint(i)->setActuatorType(dart::dynamics::Joint::ACCELERATION);
      hrp4->getJoint(i)->setForceUpperLimit(0,  forceLimit);
      hrp4->getJoint(i)->setForceLowerLimit(0, -forceLimit);
      hrp4->getJoint(i)->setPositionLimitEnforced(true);
      std::cout << hrp4->getJoint(i)->getPositionLowerLimit(0)
                << " < " << hrp4->getJoint(i)->getName()
                << " < " << hrp4->getJoint(i)->getPositionUpperLimit(0) << std::endl;
    }

  }
  // Creating the clone and its kinematics object
  Kinematics robot_kinematics(hrp4, urdf_filename);
  // initializing the controller
  auto controller = std::make_shared<Controller>(hrp4, robot_kinematics);
  

  // Wrap a WorldNode around it
  osg::ref_ptr<WorldNode> node = new WorldNode(world, hrp4, controller);
  node->setNumStepsPerCycle(5);

  // Create a Viewer and set it up with the WorldNode
  dart::gui::osg::ImGuiViewer viewer;
  viewer.addWorldNode(node);

  viewer.switchHeadlights(false);

  // Set recording
  //const std::string video_dir = "./video";
  //std::filesystem::create_directory(std::filesystem::absolute(video_dir));
  //viewer.record(video_dir,"frame");

  // Set the dimensions for the window
  viewer.setUpViewInWindow(0, 0, 1600, 900);

  // Set the window name
  viewer.realize();
  osgViewer::Viewer::Windows windows;
  viewer.getWindows(windows);
  windows.front()->setWindowName("HRP4 MPC");

  // read .txt files to store trajectory if needed - first store the file size
  std::ifstream inFile; inFile.open("./view.txt");
  if (!inFile) {
    std::cout << "Unable to open File" << std::endl;
    exit(1);
  }
  int file_size = 0; double buffer;
  while (inFile >> buffer) {
    file_size = file_size + 1;
  }
  inFile.close();
  Eigen::VectorXd viewpoint = Eigen::VectorXd::Zero(10);
  inFile.open("./view.txt");
  if (!inFile) {
    std::cout << "Unable to open File" << std::endl;
    exit(1);
  }
  for (int i = 0; i < 10; i++) {
    inFile >> viewpoint(i);
  }
  inFile.close();

  viewer.getCameraManipulator()->setHomePosition(
      ::osg::Vec3d( viewpoint(0),  viewpoint(1), viewpoint(2))*viewpoint(3),
      ::osg::Vec3d( viewpoint(4),  viewpoint(5), viewpoint(6)),
      ::osg::Vec3d( viewpoint(7),  viewpoint(8), viewpoint(9)));

  /*
  //Adjust the viewpoint of the Viewer
  viewer.getCameraManipulator()->setHomePosition(
        ::osg::Vec3d( 0.8,  -8.2*3.28*0.2, 3.3*0.155)*1.0,  
        ::osg::Vec3d( -0.10,  2.5, 0.35),  
        ::osg::Vec3d( 0.00,  0.2, 2.1));
  */
  // We need to re-dirty the CameraManipulator by passing it into the viewer
  // again, so that the viewer knows to update its HomePosition setting
  viewer.setCameraManipulator(viewer.getCameraManipulator());

  // Begin running the application loop
  viewer.simulate(true);
  try {
    viewer.run();
  }
  catch (int result) {
    std::cout << result << std::endl;
    return result;
  }
//  catch (...) {
//    std::cout << "crashed" << std::endl;
//  }

  return 1;
}
