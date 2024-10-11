#include "Kinematics.hpp"


Kinematics::Kinematics(const dart::dynamics::SkeletonPtr& robot, const std::string& urdf_filename)  {

  // List of joints to keep unlocked
  //const std::vector<std::string> joint_names = {"root_joint", "CHEST_P", "CHEST_Y", "L_SHOULDER_P", "L_SHOULDER_R", "L_SHOULDER_Y", "L_ELBOW_P", "L_WRIST_Y", "L_WRIST_P", "L_WRIST_R", "NECK_Y", "NECK_P", "R_SHOULDER_P", "R_SHOULDER_R", "R_SHOULDER_Y", "R_ELBOW_P", "R_WRIST_Y", "R_WRIST_P", "R_WRIST_R", "L_HIP_Y", "L_HIP_R", "L_HIP_P", "L_KNEE_P", "L_ANKLE_P", "L_ANKLE_R", "R_HIP_Y", "R_HIP_R", "R_HIP_P", "R_KNEE_P", "R_ANKLE_P", "R_ANKLE_R"};
  //const std::vector<std::string> joint_names = {"root_joint", "CHEST_P", "CHEST_Y", "L_SHOULDER_P", "L_SHOULDER_R", "L_SHOULDER_Y", "L_ELBOW_P", "R_SHOULDER_P", "R_SHOULDER_R", "R_SHOULDER_Y", "R_ELBOW_P", "L_HIP_Y", "L_HIP_R", "L_HIP_P", "L_KNEE_P", "L_ANKLE_P", "L_ANKLE_R", "R_HIP_Y", "R_HIP_R", "R_HIP_P", "R_KNEE_P", "R_ANKLE_P", "R_ANKLE_R"};
  //const std::vector<std::string> joint_names = {"root_joint",  "L_SHOULDER_P", "L_SHOULDER_R", "L_SHOULDER_Y", "L_ELBOW_P", "R_SHOULDER_P", "R_SHOULDER_R", "R_SHOULDER_Y", "R_ELBOW_P", "L_HIP_Y", "L_HIP_R", "L_HIP_P", "L_KNEE_P", "L_ANKLE_P", "L_ANKLE_R", "R_HIP_Y", "R_HIP_R", "R_HIP_P", "R_KNEE_P", "R_ANKLE_P", "R_ANKLE_R"};

  //  chest + minimal arms + legs
  const std::vector<std::string> joint_names = {"root_joint", "CHEST_P", "CHEST_Y", "L_SHOULDER_P", "L_SHOULDER_R", "L_ELBOW_P", "R_SHOULDER_P","R_SHOULDER_R", "R_ELBOW_P", "L_HIP_Y", "L_HIP_R", "L_HIP_P", "L_KNEE_P", "L_ANKLE_P", "L_ANKLE_R", "R_HIP_Y", "R_HIP_R", "R_HIP_P", "R_KNEE_P", "R_ANKLE_P", "R_ANKLE_R"};
  // minimal arms + legs
  //const std::vector<std::string> joint_names = {"root_joint", "L_SHOULDER_P", "L_ELBOW_P", "R_SHOULDER_P", "R_ELBOW_P", "L_HIP_Y", "L_HIP_R", "L_HIP_P", "L_KNEE_P", "L_ANKLE_P", "L_ANKLE_R", "R_HIP_Y", "R_HIP_R", "R_HIP_P", "R_KNEE_P", "R_ANKLE_P", "R_ANKLE_R"};
  // legs only
  //const std::vector<std::string> joint_names = {"root_joint", "L_HIP_Y", "L_HIP_R", "L_HIP_P", "L_KNEE_P", "L_ANKLE_P", "L_ANKLE_R", "R_HIP_Y", "R_HIP_R", "R_HIP_P", "R_KNEE_P", "R_ANKLE_P", "R_ANKLE_R"};

  // Number of joints excluding the 6-dimensional floating base.
  // Since normal joints are 1-dimensional this has same dimension
  // of the reduced configuration space minus the floating base.
  n_joints = joint_names.size() - 1;
  // Dimension of joint configuration space of the Dart model
  // including the 6-dimensional floating base.
  n_joints_full = robot->getPositions().rows();
  // Map of the joints from the list into the joint indexes of Dart.
  // Iteration starts from 1 to exclude the root_joint.
  for (int i = 1; i < joint_names.size(); ++i) {
    int index = robot->getDof(joint_names[i])->getIndexInSkeleton();
    dof_indices.push_back(index);
  }

  robot_model = buildModelFromJoints(urdf_filename, joint_names);
  current_data = pinocchio::Data(robot_model);

  r_sole_ID = robot_model.getFrameId("r_sole");
  l_sole_ID = robot_model.getFrameId("l_sole");
  // The torso will be used as the orientation to be regulated
  torso_ID = robot_model.getFrameId("torso");
  //CoM_ID = robot_model.getJointId("root_joint"); // True CoM is used instead
  r_hand_ID = robot_model.getFrameId("R_HAND_J1_LINK");
  l_hand_ID = robot_model.getFrameId("L_HAND_J1_LINK");


  J_supp2fb_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  J_swing_foot_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  J_CoM_ = Eigen::MatrixXd::Zero(3, 6+n_joints);
  J_torso_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  J_left_foot_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  J_right_foot_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  J_CMM_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  J_right_hand_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  J_left_hand_ = Eigen::MatrixXd::Zero(6, 6+n_joints);

  Jdot_supp2fb_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  Jdot_swing_foot_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  Jdot_CoM_ = Eigen::MatrixXd::Zero(3, 6+n_joints);
  Jdot_torso_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  Jdot_left_foot_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  Jdot_right_foot_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  Jdot_CMM_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  Jdot_right_hand_ = Eigen::MatrixXd::Zero(6, 6+n_joints);
  Jdot_left_hand_ = Eigen::MatrixXd::Zero(6, 6+n_joints);


  std::cout << "List of joints in the reduced model:" << std::endl;
  for (pinocchio::JointIndex joint_id = 0; joint_id < robot_model.joints.size(); ++joint_id)
    std::cout << "\t" << joint_id << " - " << robot_model.names[joint_id] << std::endl;
}


Eigen::VectorXd Kinematics::predictNextConfig() {
    Eigen::VectorXd next_config(7+n_joints);
    pinocchio::integrate(robot_model, current_config, timeStep*current_velocity, next_config);
    return next_config;
}

Eigen::VectorXd Kinematics::predictNextConfig(const Eigen::VectorXd & qdot) {
  Eigen::VectorXd next_config(7+n_joints);
  pinocchio::integrate(robot_model, current_config, timeStep*qdot, next_config);
  return next_config;
}


void Kinematics::updateKinematics(const Eigen::VectorXd & config, const Eigen::VectorXd & robot_velocity, Foot support_foot) {
  // Set the current config and compute the relevant kinematic quantities.
  // This function has to be called before getting the kinematics and must
  // be extended to compute additional quantities if required.

  current_config = config;
  current_velocity = robot_velocity;

  pinocchio::forwardKinematics(robot_model, current_data, current_config, current_velocity);

  //pinocchio::centerOfMass(robot_model, current_data, false);

  // Frame placements have to be updated separately or all together
  pinocchio::updateFramePlacement(robot_model,current_data,l_sole_ID);
  pinocchio::updateFramePlacement(robot_model,current_data,r_sole_ID);
  pinocchio::updateFramePlacement(robot_model,current_data,torso_ID);
  pinocchio::updateFramePlacement(robot_model,current_data,r_hand_ID);
  pinocchio::updateFramePlacement(robot_model,current_data,l_hand_ID);

  this->computeJacobians(support_foot);
  this->computeTaskPositions(support_foot);

}

void Kinematics::computeJacobians(Foot support_foot){

  //pinocchio::computeJointJacobians(robot_model, current_data , current_config);
  // This also computes all the jacobians apparently
  pinocchio::jacobianCenterOfMass(robot_model, current_data, false);

  J_CMM_ = pinocchio::ccrba(robot_model, current_data, current_config, current_velocity);
  pinocchio::computeCentroidalMapTimeVariation(robot_model, current_data, current_config, current_velocity);
  //pinocchio::computeJointJacobiansTimeVariation(robot_model, current_data, current_config, current_velocity);

  // Jacobians of the two feet
  pinocchio::getFrameJacobian(robot_model,current_data, r_sole_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, J_right_foot_);
  pinocchio::getFrameJacobian(robot_model, current_data, l_sole_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, J_left_foot_);
  // Jacobian of torso
  pinocchio::getFrameJacobian(robot_model, current_data, torso_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, J_torso_);
  // Jacobians of the hands
  pinocchio::getFrameJacobian(robot_model,current_data, r_hand_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, J_right_hand_);
  pinocchio::getFrameJacobian(robot_model, current_data, l_hand_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, J_left_hand_);

  J_CoM_ = current_data.Jcom;

  Jdot_CMM_ = current_data.dAg;
  Eigen::Matrix6d Ig = current_data.Ig;
  Jdot_CoM_ = (Ig.inverse()*Jdot_CMM_).topRows(3);

  pinocchio::getFrameJacobianTimeVariation(robot_model, current_data, r_sole_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, Jdot_right_foot_);
  pinocchio::getFrameJacobianTimeVariation(robot_model, current_data, l_sole_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, Jdot_left_foot_);
  pinocchio::getFrameJacobianTimeVariation(robot_model, current_data, torso_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, Jdot_torso_);
  pinocchio::getFrameJacobianTimeVariation(robot_model,current_data, r_hand_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, Jdot_right_hand_);
  pinocchio::getFrameJacobianTimeVariation(robot_model, current_data, l_hand_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED, Jdot_left_hand_);

  if(support_foot == Foot::RIGHT){
    J_swing_foot_ = J_left_foot_;
    Jdot_swing_foot_ = Jdot_left_foot_;
  }
  else{
    J_swing_foot_ = J_right_foot_;
    Jdot_swing_foot_ = Jdot_right_foot_;
  }

}

Eigen::MatrixXd Kinematics::getTaskJacobian(const Task task) const{
  // Get the Jacobians in the LOCAL_WORLD_ALIGNED frame, assuming
  // that they have been computed inside computeJacobians().
  Eigen::MatrixXd jacobian;
  switch (task) {
    case Task::COM_POSITION: {
      jacobian = J_CoM_;
      break;
    }
    case Task::TORSO_POSE: {
      jacobian = J_torso_;
      break;
    }
    case Task::TORSO_ORIENTATION: {
      jacobian = J_torso_.bottomRows(3);
      break;
    }
    case Task::SWING_FOOT_POSE: {
      jacobian = J_swing_foot_;
      break;
    }
    case Task::LEFT_FOOT_POSE: {
      jacobian = J_left_foot_;
      break;
    }
    case Task::RIGHT_FOOT_POSE: {
      jacobian = J_right_foot_;
      break;
    }
    case Task::CMM: {
      jacobian = J_CMM_;
      break;
    }
    case Task::LEFT_HAND_POSE: {
      jacobian = J_left_hand_;
      break;
    }
    case Task::RIGHT_HAND_POSE: {
      jacobian = J_right_hand_;
      break;
    }
  }
  return jacobian;
}

Eigen::MatrixXd Kinematics::getTaskJacobianTimeDeriv(const Task task) const{
  // Get the Jacobians in the LOCAL_WORLD_ALIGNED frame, assuming
  // that they have been computed inside computeJacobians().
  Eigen::MatrixXd jacobian;
  switch (task) {
    case Task::COM_POSITION: {
      jacobian = Jdot_CoM_;
      break;
    }
    case Task::TORSO_POSE: {
      jacobian = Jdot_torso_;
      break;
    }
    case Task::TORSO_ORIENTATION: {
      jacobian = Jdot_torso_.bottomRows(3);
      break;
    }
    case Task::SWING_FOOT_POSE: {
      jacobian = Jdot_swing_foot_;
      break;
    }
    case Task::LEFT_FOOT_POSE: {
      jacobian = Jdot_left_foot_;
      break;
    }
    case Task::RIGHT_FOOT_POSE: {
      jacobian = Jdot_right_foot_;
      break;
    }
    case Task::CMM: {
      jacobian = Jdot_CMM_;
      break;
    }
    case Task::LEFT_HAND_POSE: {
      jacobian = Jdot_left_hand_;
      break;
    }
    case Task::RIGHT_HAND_POSE: {
      jacobian = Jdot_right_hand_;
      break;
    }
  }
  return jacobian;
}

void Kinematics::computeTaskPositions(Foot support_foot){

  task_positions_.support_foot = support_foot;

  task_positions_.CoM.pos = current_data.com[0];
  task_positions_.CoM.ang_pos = utils::getRPY(current_data.oMf[torso_ID].rotation());

  Eigen::Vector6d torso_spatial_vel = pinocchio::getFrameVelocity(robot_model, current_data,torso_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED);
  task_positions_.CoM.ang_vel = torso_spatial_vel.tail(3);
  task_positions_.CoM.vel = current_data.vcom[0];

  // feet
  task_positions_.leftFoot.pos = current_data.oMf[l_sole_ID].translation();
  task_positions_.leftFoot.ang_pos = utils::getRPY(current_data.oMf[l_sole_ID].rotation());

  Eigen::Vector6d left_spatial_vel = pinocchio::getFrameVelocity(robot_model, current_data,l_sole_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED);
  task_positions_.leftFoot.vel = left_spatial_vel.head(3);
  task_positions_.leftFoot.ang_vel = left_spatial_vel.tail(3);

  task_positions_.rightFoot.pos = current_data.oMf[r_sole_ID].translation();
  task_positions_.rightFoot.ang_pos = utils::getRPY(current_data.oMf[r_sole_ID].rotation());

  Eigen::Vector6d right_spatial_vel = pinocchio::getFrameVelocity(robot_model, current_data,r_sole_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED);
  task_positions_.rightFoot.vel = right_spatial_vel.head(3);
  task_positions_.rightFoot.ang_vel = right_spatial_vel.head(3);

  // hands
  task_positions_.leftHand.pos = current_data.oMf[l_hand_ID].translation();
  task_positions_.leftHand.ang_pos = utils::getRPY(current_data.oMf[l_hand_ID].rotation());

  left_spatial_vel = pinocchio::getFrameVelocity(robot_model, current_data,l_hand_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED);
  task_positions_.leftHand.vel = left_spatial_vel.head(3);
  task_positions_.leftHand.ang_vel = left_spatial_vel.tail(3);

  task_positions_.rightHand.pos = current_data.oMf[r_hand_ID].translation();
  task_positions_.rightHand.ang_pos = utils::getRPY(current_data.oMf[r_hand_ID].rotation());

  right_spatial_vel = pinocchio::getFrameVelocity(robot_model, current_data,r_hand_ID, pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED);
  task_positions_.rightHand.vel = right_spatial_vel.head(3);
  task_positions_.rightHand.ang_vel = right_spatial_vel.head(3);
}

Eigen::VectorXd Kinematics::getTaskPosition(const Task task) const{
  Eigen::VectorXd task_pos;
  switch (task) {
    case Task::COM_POSITION: {
      task_pos = task_positions_.CoM.pos;
      break;
    }
    case Task::TORSO_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.CoM.pos, task_positions_.CoM.ang_pos;
      task_pos = pose;
      break;
    }
    case Task::TORSO_ORIENTATION: {
      task_pos = task_positions_.CoM.ang_pos;
      break;
    }
    case Task::SWING_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      if (task_positions_.support_foot == Foot::RIGHT) // Right is support, left is swing
        pose << task_positions_.leftFoot.pos, task_positions_.leftFoot.ang_pos;
      else
        pose << task_positions_.rightFoot.pos, task_positions_.rightFoot.ang_pos;
      task_pos = pose;
      break;
    }
    case Task::SUPPORT_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      if (task_positions_.support_foot == Foot::RIGHT) // Right is support, left is swing
        pose << task_positions_.rightFoot.pos, task_positions_.rightFoot.ang_pos;
      else
        pose << task_positions_.leftFoot.pos, task_positions_.leftFoot.ang_pos;
      task_pos = pose;
      break;
    }
    case Task::LEFT_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.leftFoot.pos, task_positions_.leftFoot.ang_pos;
      task_pos = pose;
      break;
    }
    case Task::RIGHT_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.rightFoot.pos, task_positions_.rightFoot.ang_pos;
      task_pos = pose;
      break;
    }
    case Task::LEFT_HAND_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.leftHand.pos, task_positions_.leftHand.ang_pos;
      task_pos = pose;
      break;
    }
    case Task::RIGHT_HAND_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.rightHand.pos, task_positions_.rightHand.ang_pos;
      task_pos = pose;
      break;
    }
    default: {
      task_pos = Eigen::Vector3d::Zero();
      break;
    }
  }
  return task_pos;
}


Eigen::VectorXd Kinematics::getTaskVelocity(const Task task) const{
  Eigen::VectorXd task_vel;
  switch (task) {
    case Task::COM_POSITION: {
      task_vel = task_positions_.CoM.vel;
      break;
    }
    case Task::TORSO_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.CoM.vel, task_positions_.CoM.ang_vel;
      task_vel = pose;
      break;
    }
    case Task::TORSO_ORIENTATION: {
      task_vel = task_positions_.CoM.ang_vel;
      break;
    }
    case Task::SWING_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      if (task_positions_.support_foot == Foot::RIGHT) // Right is support, left is swing
        pose << task_positions_.leftFoot.vel, task_positions_.leftFoot.ang_vel;
      else
        pose << task_positions_.rightFoot.vel, task_positions_.rightFoot.ang_vel;
      task_vel = pose;
      break;
    }
    case Task::SUPPORT_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      if (task_positions_.support_foot == Foot::RIGHT) // Right is support, left is swing
        pose << task_positions_.rightFoot.vel, task_positions_.rightFoot.ang_vel;
      else
        pose << task_positions_.leftFoot.vel, task_positions_.leftFoot.ang_vel;
      task_vel = pose;
      break;
    }
    case Task::LEFT_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.leftFoot.vel, task_positions_.leftFoot.ang_vel;
      task_vel = pose;
      break;
    }
    case Task::RIGHT_FOOT_POSE: {
      Eigen::VectorXd pose(6);
      pose << task_positions_.rightFoot.vel, task_positions_.rightFoot.ang_vel;
      task_vel = pose;
      break;
    }
    default: {
      task_vel = Eigen::Vector3d::Zero();
      break;
    }
  }
  return task_vel;
}

Eigen::VectorXd Kinematics::coords2joints(const Eigen::VectorXd & c){
  // Mapping between the Dart configuration to the joints
  Eigen::VectorXd q(n_joints+7);
  for(int i = 0; i < n_joints; i++){
    int index = dof_indices[i];
    q(i+7) = c(index);
  }
  q.head(7) = Kinematics::D2PFloatingBase(c.head(6));
  return q;
}

/*Eigen::MatrixXd Kinematics::coords2joints(const Eigen::MatrixXd & A){
  Eigen::MatrixXd B(n_joints, A.cols());
  for(int i = 0; i < B.cols(); i++){
    Eigen::VectorXd temp = A.col(i);
    B.col(i) = Kinematics::coords2joints(temp);
  }
  return B;
}*/

Eigen::VectorXd Kinematics::joints2coords(const Eigen::VectorXd & q){
  // Mapping between the n_joints joints to the n_joints_full coordinates
  // but it does not populate the floating base coordinates (first 6 for DART)
  Eigen::VectorXd c(n_joints_full);
  c.setZero();
  for(int i = 0; i < n_joints; i++){
    int index = dof_indices[i];
    c(index) = q(i+6); // add 6 to skip the floating base
  }
  return c;
}

Eigen::VectorXd Kinematics::D2PFloatingBase(const Eigen::Vector6d & dart_pose) {
  Eigen::Isometry3d iso = dart::dynamics::FreeJoint::convertToTransform(dart_pose);
  Eigen::Quaterniond quat(iso.rotation()); quat.normalize();

  Eigen::VectorXd pinocchio_pose = Eigen::VectorXd::Zero(7);
  pinocchio_pose.head(3) = iso.translation();
  pinocchio_pose.segment(3,4) = quat.coeffs();
  return pinocchio_pose;
}


int Kinematics::getDofs() const {return this->n_joints;}
int Kinematics::getDofsFullModel() const {return this->n_joints_full;}


void Kinematics::debug() {

}



pinocchio::Model Kinematics::buildModelFromJoints(const std::string & urdf_filename, const std::vector<std::string> & joint_names) {

  pinocchio::Model full_model;
  pinocchio::urdf::buildModel(urdf_filename,pinocchio::JointModelFreeFlyer(),full_model);


  std::vector<pinocchio::JointIndex> list_of_joints_to_keep_unlocked_by_id;
  for(auto it = joint_names.begin(); it != joint_names.end(); ++it) {
    const std::string & joint_name = *it;
    if(full_model.existJointName(joint_name))
      list_of_joints_to_keep_unlocked_by_id.push_back(full_model.getJointId(joint_name));
    else std::cout << "joint: " << joint_name << " does not belong to the model";
  }
  // Transform the list into a list of joints to lock
  std::vector<pinocchio::JointIndex> list_of_joints_to_lock_by_id;
  for(pinocchio::JointIndex joint_id = 1; joint_id < full_model.joints.size(); ++joint_id) {
    const std::string joint_name = full_model.names[joint_id];
    if(utils::is_in_vector(joint_names,joint_name))
      continue;
    else list_of_joints_to_lock_by_id.push_back(joint_id);
  }

  Eigen::VectorXd q_full = D2PInitialConfig(full_model);
  // Build the reduced model from the list of lock joints
  return buildReducedModel(full_model,list_of_joints_to_lock_by_id,q_full);
}


Eigen::VectorXd Kinematics::D2PInitialConfig(const pinocchio::Model & model) { // TODO joint values should be written automatically



  Eigen::VectorXd config = pinocchio::neutral(model);

  // Eigen::Vector6d free_positions; // Floating base pose using DART's coordinates
  // free_positions << 0.0 , 4*M_PI/180.0, 0.0 , 0.035502257,-1.5, 0.751;
  //config.head(7) = D2PFloatingBase(free_positions);

  // Note that we have to shift by 5 to account for the floating base part but discounting the universe (0) and root_joint (1)
  // right leg
  config(model.getJointId("R_HIP_Y")+5)= 0.0;
  config(model.getJointId("R_HIP_R")+5)=  -3*M_PI/180 ;
  config(model.getJointId("R_HIP_P")+5)=  -25*M_PI/180 ;
  config(model.getJointId("R_KNEE_P")+5)= 50*M_PI/180 ;
  config(model.getJointId("R_ANKLE_P")+5)=  -30*M_PI/180 ;
  config(model.getJointId("R_ANKLE_R")+5)=  4*M_PI/180 ;
  // left leg
  config(model.getJointId("L_HIP_Y")+5)=  0.0 ;
  config(model.getJointId("L_HIP_R")+5)=  3*M_PI/180 ;
  config(model.getJointId("L_HIP_P")+5)=  -25*M_PI/180 ;
  config(model.getJointId("L_KNEE_P")+5)=  50*M_PI/180 ;
  config(model.getJointId("L_ANKLE_P")+5)=  -30*M_PI/180 ;
  config(model.getJointId("L_ANKLE_R")+5)=  -4*M_PI/180 ;
  // right arm
  config(model.getJointId("R_SHOULDER_P")+5)=  (4)*M_PI/180 ;
  config(model.getJointId("R_SHOULDER_R")+5)=  -8*M_PI/180  ;
  config(model.getJointId("R_SHOULDER_Y")+5)=  0 ;
  config(model.getJointId("R_ELBOW_P")+5)=  -25*M_PI/180 ;
  // left arm
  config(model.getJointId("L_SHOULDER_P")+5)=  (4)*M_PI/180 ;
  config(model.getJointId("L_SHOULDER_R")+5)=  8*M_PI/180 ;
  config(model.getJointId("L_SHOULDER_Y")+5)=  0 ;
  config(model.getJointId("L_ELBOW_P")+5)=  -25*M_PI/180 ;

  return config;
}
