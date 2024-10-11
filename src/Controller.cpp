#include "Controller.hpp"

#include <utility>
#include <chrono>

Controller::Controller(dart::dynamics::SkeletonPtr _robot, Kinematics _robot_kinematics)
    : mRobot(std::move(_robot)), robotKinematics(std::move(_robot_kinematics))
{
  std::cout << "[Controller] initialization" << std::endl;

  this->setInitialConfiguration();

  this->setParamsFromFile("./simulation_params.toml");

  mLeftFoot = mRobot->getBodyNode("l_sole");
  mRightFoot = mRobot->getBodyNode("r_sole");
  mTorso = mRobot->getBodyNode("torso");

  // set the zmp at the center of the support polygon
  zmp_buffer_box = 0.5*(mRightFoot->getWorldTransform().translation() + mLeftFoot->getWorldTransform().translation());

//  std::cout << " QUANTITIES FROM DART:"  << std::endl;
//  std::cout << "Total mass = " << mRobot->getMass() << std::endl;
//  std::cout << "CoM position = " << mRobot->getCOM().transpose() << std::endl;
//  std::cout << "LEFT foot = " << mLeftFoot->getCOM().transpose() << std::endl;
//  std::cout << "RIGHT foot = " << mRightFoot->getCOM().transpose() << std::endl;

  walk_state_.global_iter = 0;
  walk_state_.mpc_iter = 0;
  walk_state_.footstep_idx = 0;
  /* Set the first support foot, but since there is a dummy footstep to move the
   * ZMP from the center of the double support polygon to the next support foot
   * this will also be the first actual swing foot. */
  walk_state_.support_foot = Foot::RIGHT;
  walk_state_.mode = WalkMode::STANDING;

  Foot pred_support_foot = walk_state_.support_foot;
  for(int i = 0; i<N+1; ++i) {
    if(i % (S+D) == 0 && i != 0) switchSupportFoot(pred_support_foot);
    walk_state_.predicted_support_foot.push_back(pred_support_foot);
  }

  // Initialize desired state with reasonable values
  desired.CoM.pos = Eigen::Vector3d(0.0, 0.0, comTargetHeight);
  desired.CoM.vel = Eigen::Vector3d::Zero();
  desired.CoM.ang_pos = Eigen::Vector3d::Zero();//
  desired.CoM.ang_vel = Eigen::Vector3d::Zero();


  n_joints = robotKinematics.getDofs();
  n_joints_full = robotKinematics.getDofsFullModel();


  predictedJointPos.resize(7+n_joints,N+1);
  predictedCoMPos.resize(3,N+1);
  swing_pin.resize(6,N);

  QAcc = Eigen::MatrixXd::Zero(6+n_joints, N+1);


  // Fill matrix of joint limits for all the joints, then
  Eigen::VectorXd joint_lims_full_lb = Eigen::VectorXd::Zero(n_joints_full);
  Eigen::VectorXd joint_lims_full_ub = Eigen::VectorXd::Zero(n_joints_full);

  for (int i = 0; i < mRobot->getNumJoints(); i++) {
    size_t dim = mRobot->getJoint(i)->getNumDofs();
    if(dim==1) {
      int joint_idx = mRobot->getDof(mRobot->getJoint(i)->getName())->getIndexInSkeleton();
      joint_lims_full_lb(joint_idx) = mRobot->getJoint(i)->getPositionLowerLimit(0);
      joint_lims_full_ub(joint_idx) = mRobot->getJoint(i)->getPositionUpperLimit(0);
    }
  }
  Eigen::MatrixXd joint_lims = Eigen::MatrixXd::Zero(n_joints, 2);
  joint_lims.col(0) = robotKinematics.coords2joints(joint_lims_full_lb).tail(n_joints);
  joint_lims.col(1) = robotKinematics.coords2joints(joint_lims_full_ub).tail(n_joints);


  mpc = std::make_unique<JointLevelISMPC>(n_joints, joint_lims, getCurrentConfiguration());

  //FIXME fill initial state for ismpc, with support foot bool
  initial.CoM.pos = mRobot->getCOM();
  initial.CoM.vel = Eigen::Vector3d::Zero();
  initial.zmpPos = mRobot->getCOM();
  initial.zmpPos(2) = 0.0;

  initial.leftFoot.pos = mLeftFoot->getCOM();
  initial.leftFoot.ang_pos = Eigen::Vector3d::Zero();
  initial.rightFoot.pos = mRightFoot->getCOM();
  initial.rightFoot.ang_pos = Eigen::Vector3d::Zero();

  filteredZmpPos = initial.CoM.pos;


  mpc_data_.current_state = initial;

  planFootsteps();


  // Can be used to compute a rotated the initial pose in dart coordinates
  //Eigen::Vector6d rot_pose = dart::dynamics::FreeJoint::convertToPositions(rotz(M_PI_2)*dart::dynamics::FreeJoint::convertToTransform(mRobot->getPositions().head(6)));

  this->setLogger();
}

Controller::~Controller() = default;

void Controller::update() {

  if (walk_state_.global_iter == (int)((num_steps+1)*(singleSupportDuration+doubleSupportDuration)/worldTimeStep)) {
    std::cout << "exiting" << std::endl;
    if(mRobot->getCOM()[2] < 0.45) {
      throw 1;
    }
    else {
      throw 0;
    }

  }

//  std::cout << "[Controller] update" << std::endl;

  if (walk_state_.global_iter>=push_timing_ && walk_state_.global_iter<= push_timing_+push_duration_) 
    mTorso->addExtForce(push_magnitude_*Eigen::Vector3d(-cos(push_angle_/180.0*M_PI),-sin(push_angle_/180.0*M_PI),0.0));

  // Draw on the screen
  draw();

  ComAccXYMeas = (mRobot->getCOMLinearAcceleration()).head(2);
  // Retrieve current state
  current.CoM.pos = mRobot->getCOM();
  current.CoM.vel = mRobot->getCOMLinearVelocity();
  current.zmpPos = current.CoM.pos - mRobot->getCOMLinearAcceleration() / (omega*omega);
  current.zmpPos = getZmpFromExternalForces();
  jointPos = mRobot->getPositions();
  jointVel = mRobot->getVelocities();
  current.leftFoot.pos = mLeftFoot->getCOM();
  current.leftFoot.ang_pos = utils::getRPY(mLeftFoot);
  current.rightFoot.pos = mRightFoot->getCOM();
  current.rightFoot.ang_pos = utils::getRPY(mRightFoot);

  // mpc_data_.current_state = current;
  mpc_data_.current_q = robotKinematics.coords2joints(jointPos);
  mpc_data_.current_qdot = this->getCurrentVelocity();

  // ZMP measure
  filteredZmpPos = getZmpFromExternalForces();
  mpc_data_.filtered_zmp_pos = filteredZmpPos;

  // Timer
  auto start = std::chrono::high_resolution_clock::now();

  if (walk_state_.global_iter == 0) {
    walk_state_.mode = WalkMode::WALKING;
  }

  // Compute predictions and references
  this->computePredictions();
  this->computeReferences();



  // Prepare MPC matrices
  mpc->prepareMPC(mpc_data_, walk_state_);

  // Feasibility check
  auto [bxl, bxu, byl, byu] = mpc->getFeasibilityRegionBounds(mpc_data_);
  feasibility_bounds << bxl, bxu, byl, byu; // TODO see if this is needed

//  double epsilon = 0.1;
//  Eigen::Vector3d pk = mpc_data_.CoM_preds.pos[0] + mpc_data_.CoM_preds.jacobians[0] * mpc_data_.qdot[0] * (1/omega + predictionTime);
//  bool is_feasible = (pk(0) < bxu - epsilon) && (pk(0) > bxl + epsilon) && (pk(1) < byu - epsilon) && (pk(1) > byl + epsilon);

  /*if(walk_state_.global_iter==395) { //!is_feasible) {
    std::cout << "HRP4 is about to die" << std::endl;
    for(int i = walk_state_.footstep_idx + 2; i < (footstep_planner->footstepPlan).size(); ++i) {
      footstep_planner->footstepPlan.at(i)(0) += 0.1;
    }
    footstep_planner->computeModeSequence();
    footstep_planner->computeSwingTrajectory();
    footstep_planner->computeZMPTrajectory();
    this->computeReferences();
    mpc->prepareMPC(mpc_data_, walk_state_);
    auto [bxl, bxu, byl, byu] = mpc->getFeasibilityRegionBounds(mpc_data_);
    feasibility_bounds = Eigen::Vector4d(bxl, bxu, byl, byu);
    std::cout << "re-bounds: " << feasibility_bounds.transpose() << std::endl;
  }*/

  swing_des = footstep_planner->getSwingPos(walk_state_.global_iter);
  
  // Solve MPC QP
  QAcc = mpc->computeQAcc(mpc_data_, walk_state_);

  Eigen::VectorXd qAcc_coords = robotKinematics.joints2coords(QAcc.col(0));
  for(int i = 6; i < qAcc_coords.rows(); i++) mRobot->setCommand(i, qAcc_coords(i));


  // UPDATE INDEXES
  ++walk_state_.global_iter;

  if(walk_state_.global_iter % (singleSupportSamples+doubleSupportSamples) == 0 && walk_state_.global_iter !=0) {
    switchSupportFoot(walk_state_.support_foot);
    std::cout << "\033[1;32m[Controller] set support foot to "
              << supportFoot2String(walk_state_.support_foot)
              << " at time "
              << walk_state_.global_iter*worldTimeStep
              << " s \033[0m"  <<std::endl;
  }

  if(walk_state_.global_iter % static_cast<int>(timeStep/worldTimeStep) == 0) {

    ++walk_state_.mpc_iter;
    if(walk_state_.mpc_iter == (S+D)) {
      walk_state_.mpc_iter = 0;
    }
    predictLastSupportFoot();
  }

  walk_state_.footstep_idx = floor(walk_state_.global_iter*worldTimeStep/(singleSupportDuration+doubleSupportDuration));


  // Shift solution for next step
  //QDots.leftCols(N-1) = QDots_.rightCols(N-1);
  //QDots.rightCols(1) = QDots_.rightCols(1); // repeats last command

  // TIMINGS
  // auto finish = std::chrono::high_resolution_clock::now();
  // auto interval = std::chrono::time_point_cast<std::chrono::microseconds>(finish) - std::chrono::time_point_cast<std::chrono::microseconds>(start);
  // std::cout << interval.count()/1000.0 << " ms" << std::endl;

  //LOGGING
  // Store the results in files (for plotting)

  // Fill prediction matrices for logging
  /*for(int j = 0; j <= N; j++) {
    predictedJointPos.col(j) = mpc_data_.q[j];
    //predictedCoMPos.col(j) = mpc_data_.CoM_preds.pos[j].head(3);
  }*/
  //com_jacobian = mpc_data_.CoM_preds.jacobians[0].topRows(2).transpose();
  right_hand_pos_log = mpc_data_.right_hand_preds.pos[0];
  right_hand_pos_ref_log = mpc_data_.right_hand_refs.pos[0];

  ZmpRef = mpc_data_.zmp_refs.pos[0];

  logger_.logList();

  mpc_data_.clearAll();

}

void Controller::planFootsteps() {

  int total_samples = num_steps*(singleSupportSamples+doubleSupportSamples);
  std::vector<Vref> vrefs(total_samples, vref_const);
  for(int i = 0; i<(singleSupportSamples+doubleSupportSamples);++i) vrefs[i] = Vref(0.0,0.0,0.0);
  for(int i = total_samples-3*(singleSupportSamples+doubleSupportSamples);
       i<total_samples-2*(singleSupportSamples+doubleSupportSamples);++i) vrefs[i] = Vref(0.5*vref_const.x,0.5*vref_const.y,0.5*vref_const.omega);
  for(int i = total_samples-2*(singleSupportSamples+doubleSupportSamples);
      i<total_samples;++i) vrefs[i] = Vref(0.0,0.0,0.0);
  
  //std::vector<Vref> vrefs;
  //for(int i = 0; i<3*(singleSupportSamples+doubleSupportSamples);++i) vrefs.push_back(Vref(0.0,0.0,0.0));
  //for(int i = 0; i<2*(singleSupportSamples+doubleSupportSamples);++i) vrefs.push_back(vref_const);
  //for(int i = 0; i<2*(singleSupportSamples+doubleSupportSamples);++i) vrefs.push_back(Vref(0.0,0.0,0.0));
  //for(int i = 0; i<50*(singleSupportSamples+doubleSupportSamples);++i) vrefs.push_back(vref_const);
  
  //for(int i = 0; i<50*(singleSupportSamples+doubleSupportSamples);++i) vrefs.push_back(Vref(0.0,0.0,0.0));

  Eigen::Vector6d start_left; start_left << mLeftFoot->getCOM(), utils::getRPY(mLeftFoot);
  Eigen::Vector6d start_right; start_right << mRightFoot->getCOM(), utils::getRPY(mRightFoot);

  footstep_planner = std::make_shared<FootstepPlanner>();
  footstep_planner->plan(vrefs, start_left, start_right, walk_state_.support_foot, walk_state_.global_iter);
  std::cout << "Planned " <<footstep_planner->getSize() << " steps at positions:" << std::endl;
  for(auto item : footstep_planner->getPlan()) std::cout << item.transpose() << std::endl;

  footstep_planner->computeSwingTrajectory();
  footstep_planner->computeZMPTrajectory();
  
  auto plan_file = std::ofstream(realpath("./data/", nullptr) + std::string("/FootstepPlan.csv"), std::ofstream::out);
  for(auto item : footstep_planner->getPlan()) plan_file << item.transpose().format(utils::CSVFormat) << std::endl;
  plan_file.close();

}

void Controller::computePredictions(){

  Eigen::VectorXd q = this->getCurrentConfiguration();
  Eigen::VectorXd qDot(n_joints+6);
  qDot = this->getCurrentVelocity();

  for(int i = 0; i < N; ++i){

    robotKinematics.updateKinematics(q, qDot, walk_state_.predicted_support_foot.at(i));

    mpc_data_.CoM_preds.jacobians.push_back(
        robotKinematics.getTaskJacobian(Task::COM_POSITION));
    mpc_data_.torso_orientation_preds.jacobians.push_back(
        robotKinematics.getTaskJacobian(Task::TORSO_ORIENTATION));
    mpc_data_.left_foot_preds.jacobians.push_back(
        robotKinematics.getTaskJacobian(Task::LEFT_FOOT_POSE));
    mpc_data_.right_foot_preds.jacobians.push_back(
        robotKinematics.getTaskJacobian(Task::RIGHT_FOOT_POSE));
    mpc_data_.right_hand_preds.jacobians.emplace_back(
        robotKinematics.getTaskJacobian(Task::RIGHT_HAND_POSE));

    mpc_data_.CoM_preds.jacobians_dot.push_back(
        robotKinematics.getTaskJacobianTimeDeriv(Task::COM_POSITION));
    mpc_data_.torso_orientation_preds.jacobians_dot.push_back(
        robotKinematics.getTaskJacobianTimeDeriv(Task::TORSO_ORIENTATION));
    mpc_data_.left_foot_preds.jacobians_dot.push_back(
        robotKinematics.getTaskJacobianTimeDeriv(Task::LEFT_FOOT_POSE));
    mpc_data_.right_foot_preds.jacobians_dot.push_back(
        robotKinematics.getTaskJacobianTimeDeriv(Task::RIGHT_FOOT_POSE));
    mpc_data_.right_hand_preds.jacobians_dot.emplace_back(
        robotKinematics.getTaskJacobianTimeDeriv(Task::RIGHT_HAND_POSE));

    if(i==0){

      ComAccXY = (mpc_data_.CoM_preds.jacobians[0]*QAcc.col(0) + mpc_data_.CoM_preds.jacobians_dot[0]*qDot).head(2);
      State filtered = updateKF(mpc_data_.current_state, current, ComAccXY);
//      State filtered = updateKF(mpc_data_.current_state, current, ComAccXYMeas);
      mpc_data_.current_state = filtered;

      mpc_data_.current_state.CoM.ang_vel = robotKinematics.getTaskJacobian(Task::CMM).bottomRows(3)*qDot; // This is actually the angular momentum!


//      mpc_data_.CoM_preds.pos.emplace_back(
//          robotKinematics.getTaskPosition(Task::COM_POSITION));
      mpc_data_.torso_orientation_preds.pos.emplace_back(
          robotKinematics.getTaskPosition(Task::TORSO_ORIENTATION));
      mpc_data_.left_foot_preds.pos.emplace_back(
          robotKinematics.getTaskPosition(Task::LEFT_FOOT_POSE));
      mpc_data_.right_foot_preds.pos.emplace_back(
          robotKinematics.getTaskPosition(Task::RIGHT_FOOT_POSE));
      mpc_data_.right_hand_preds.pos.emplace_back(
          robotKinematics.getTaskPosition(Task::RIGHT_HAND_POSE));
    }

    mpc_data_.centroidal_momentum_preds.jacobians.push_back(robotKinematics.getTaskJacobian(Task::CMM).bottomRows(3));
    mpc_data_.centroidal_momentum_preds.jacobians_dot.push_back(robotKinematics.getTaskJacobianTimeDeriv(Task::CMM).bottomRows(3));

    mpc_data_.q.push_back(q);
    mpc_data_.qdot.push_back(qDot);

    q = robotKinematics.predictNextConfig(qDot + timeStep*QAcc.col(i)/2);
    //q = robotKinematics.predictNextConfig();
    qDot += timeStep*QAcc.col(i);

  }

  robotKinematics.updateKinematics(q, qDot, walk_state_.predicted_support_foot.at(N));

//  mpc_data_.CoM_preds.pos.emplace_back(
//      robotKinematics.getTaskPosition(Task::COM_POSITION));
//  mpc_data_.torso_orientation_preds.pos.emplace_back(
//      robotKinematics.getTaskPosition(Task::TORSO_ORIENTATION));
//  mpc_data_.left_foot_preds.pos.emplace_back(
//      robotKinematics.getTaskPosition(Task::LEFT_FOOT_POSE));
//  mpc_data_.right_foot_preds.pos.emplace_back(
//      robotKinematics.getTaskPosition(Task::RIGHT_FOOT_POSE));
//  mpc_data_.right_hand_preds.pos.emplace_back(
//      robotKinematics.getTaskPosition(Task::RIGHT_HAND_POSE));

  mpc_data_.q.push_back(q);

}


void Controller::computeReferences(){
  mpc_data_.clearReferences();

  int predicted_footstep_idx = walk_state_.footstep_idx;

  for(int i = 0; i < N; i++){

    if ((walk_state_.mpc_iter+i+1)%(S+D)==0) {
      predicted_footstep_idx++;
    }

    int t_kpi = walk_state_.global_iter + i * timesRatio;

    Eigen::Vector6d torso_pose, torso_vel;
    torso_pose.head(5) = desired.getCoMPose().head(5);
    torso_pose(5) = footstep_planner->getTorsoOrientation(t_kpi+1*timesRatio);
    torso_vel.head(5) = desired.getCoMVelocity().head(5);
    torso_vel(5) = footstep_planner->getTorsoOrientationVel(t_kpi+1*timesRatio);

    mpc_data_.CoM_refs.pos.emplace_back(torso_pose.head(3));
    mpc_data_.CoM_refs.vel.emplace_back(torso_vel.head(3));

    mpc_data_.torso_orientation_refs.pos.emplace_back(torso_pose.tail(3));
    mpc_data_.torso_orientation_refs.vel.emplace_back(torso_vel.tail(3));


    if (walk_state_.predicted_support_foot.at(i+1) == Foot::RIGHT) {
      mpc_data_.left_foot_refs.pos.emplace_back(
          footstep_planner->getSwingPos(t_kpi+1*timesRatio));
      mpc_data_.left_foot_refs.vel.emplace_back(
          footstep_planner->getSwingVel(t_kpi+1*timesRatio));
      mpc_data_.right_foot_refs.pos.emplace_back(
          footstep_planner->getFootstepPose(predicted_footstep_idx));
      mpc_data_.right_foot_refs.vel.emplace_back(Eigen::Vector6d::Zero());
    }
    else
    {
      mpc_data_.left_foot_refs.pos.emplace_back(
          footstep_planner->getFootstepPose(predicted_footstep_idx));
      mpc_data_.left_foot_refs.vel.emplace_back(Eigen::Vector6d::Zero());
      mpc_data_.right_foot_refs.pos.emplace_back(
          footstep_planner->getSwingPos(t_kpi+1*timesRatio));
      mpc_data_.right_foot_refs.vel.emplace_back(
          footstep_planner->getSwingVel(t_kpi+1*timesRatio));
    }

    Eigen::Vector6d hand_ref;
    //hand_ref << mpc_data_.CoM_preds.pos[i] + Eigen::Vector3d(0.3, -0.4, 0.2*sin(t_kpi*worldTimeStep)), Eigen::Vector3d::Zero();
    hand_ref << Eigen::Vector3d(std::max(0.0, -3*vref_const.x + vref_const.x*t_kpi*worldTimeStep),
                                -0.4 + 0.1*sin(2*t_kpi*worldTimeStep),
                                comTargetHeight + 0.0),
                                Eigen::Vector3d(vref_const.x,
                                0.0,
                                0.0);
    /*hand_ref << Eigen::Vector3d(0.2,
                                -0.2,
                                comTargetHeight + 0.1*sin(2*t_kpi*worldTimeStep)),
                                Eigen::Vector3d(vref_const.x,
                                0.0,
                                0.0);*/
    mpc_data_.right_hand_refs.pos.emplace_back(hand_ref);
    mpc_data_.right_hand_refs.vel.emplace_back(Eigen::Vector6d::Zero());


    mpc_data_.zmp_refs.pos.emplace_back(footstep_planner->getZMPPos(t_kpi));
 //   mpc_data_.zmp_refs.vel.emplace_back(footstep_planner->getZMPVel(t_kpi));
/*          mpc_data_.zmp_refs.pos.emplace_back(0.5*(mRightFoot->getWorldTransform().translation() + mLeftFoot->getWorldTransform().translation()).head(2));
        mpc_data_.zmp_refs.vel.emplace_back(Eigen::Vector2d::Zero());*/
  }

  for(int i=N; i<N+(prev+N); ++i) {
    int t_kpi = walk_state_.global_iter + i * timesRatio;
    mpc_data_.zmp_prev.emplace_back(footstep_planner->getZMPPos(t_kpi));
    //mpc_data_.zmp_prev.emplace_back(0.5*(mRightFoot->getWorldTransform().translation() + mLeftFoot->getWorldTransform().translation()).head(2));
  }

}


void Controller::predictLastSupportFoot() {
  Foot pred_support_foot = walk_state_.predicted_support_foot.at(N-(S+D)+1);
  switchSupportFoot(pred_support_foot);
  walk_state_.predicted_support_foot.push_back(pred_support_foot);
  walk_state_.predicted_support_foot.pop_front();
}

Eigen::VectorXd Controller::getCurrentConfiguration(){
  // Converts between the full DART configuration and the
  // pinocchio compatible one with reduced number of joints.
  Eigen::VectorXd c = mRobot->getPositions();
  Eigen::VectorXd q = Eigen::VectorXd::Zero(7+n_joints);
  q.head(7) = Kinematics::D2PFloatingBase(c.head(6));
  q.tail(n_joints) = robotKinematics.coords2joints(c).tail(n_joints);
  return q;
}

Eigen::VectorXd Controller::getCurrentVelocity(){
  // Converts between the full DART joint velocities and the
  // pinocchio compatible ones with reduced number of joints.
  Eigen::VectorXd cdot = mRobot->getVelocities();
  Eigen::VectorXd qdot = Eigen::VectorXd::Zero(6+n_joints);
//  qdot.head(6) = cdot.head(6);
  qdot.head(3) = cdot.segment(3,3);
  qdot.segment(3,3) = cdot.head(3);
  qdot.tail(n_joints) = robotKinematics.coords2joints(cdot).tail(n_joints);

  return qdot;
}

State Controller::updateKF(State filtered, State current, Eigen::Vector2d comAcc)
{
    // CT matrices
    Eigen::MatrixXd A_lip = Eigen::MatrixXd::Zero(2,2);
    Eigen::VectorXd B_lip = Eigen::VectorXd::Zero(2);
    A_lip << 1,worldTimeStep,0,1;
    B_lip << worldTimeStep*worldTimeStep/2,worldTimeStep;
    
    Eigen::Vector2d x_measure = Eigen::Vector2d(current.CoM.pos(0), current.CoM.vel(0));
    Eigen::Vector2d y_measure = Eigen::Vector2d(current.CoM.pos(1), current.CoM.vel(1));
    Eigen::Vector2d x_est = Eigen::Vector2d(filtered.CoM.pos(0), filtered.CoM.vel(0));
    Eigen::Vector2d y_est = Eigen::Vector2d(filtered.CoM.pos(1), filtered.CoM.vel(1));
    
    Eigen::MatrixXd F_kf = A_lip;
    Eigen::MatrixXd G_kf = B_lip;
    Eigen::MatrixXd H_kf = Eigen::Matrix2d::Identity();
    
    Eigen::MatrixXd R_kf = Eigen::MatrixXd::Identity(2,2);
    R_kf.diagonal() << cov_meas_pos, cov_meas_vel;
    Eigen::MatrixXd Q_kf = Eigen::MatrixXd::Identity(2,2);
    Q_kf.diagonal() << cov_mod_pos, cov_mod_vel;
    
    double input_x = comAcc(0);
    double input_y = comAcc(1);
    
    Eigen::VectorXd x_pred = F_kf * x_est + G_kf * input_x;
    Eigen::MatrixXd cov_x_pred = F_kf * cov_x * F_kf.transpose() + Q_kf;

    Eigen::MatrixXd K_kf = cov_x_pred * H_kf.transpose() * (H_kf * cov_x_pred * H_kf.transpose() + R_kf).inverse();

    x_est = x_pred + K_kf * (x_measure - H_kf * x_pred);
    cov_x = (Eigen::MatrixXd::Identity(2,2) - K_kf * H_kf) * cov_x_pred * (Eigen::MatrixXd::Identity(2,2) - K_kf * H_kf).transpose() + K_kf * R_kf * K_kf.transpose();

    Eigen::VectorXd y_pred = F_kf * y_est + G_kf * input_y;
    Eigen::MatrixXd cov_y_pred = F_kf * cov_y * F_kf.transpose() + Q_kf;
    
    K_kf = cov_y_pred * H_kf.transpose() * (H_kf * cov_y_pred * H_kf.transpose() + R_kf).inverse();
    
    y_est = y_pred + K_kf * (y_measure - H_kf * y_pred);
    cov_y = (Eigen::MatrixXd::Identity(2,2) - K_kf * H_kf) * cov_y_pred * (Eigen::MatrixXd::Identity(2,2) - K_kf * H_kf).transpose() + K_kf * R_kf * K_kf.transpose();
    
    filtered.CoM.pos = Eigen::Vector3d(x_est(0), y_est(0), current.CoM.pos(2));
    filtered.CoM.vel = Eigen::Vector3d(x_est(1), y_est(1), current.CoM.vel(2));
    
    return filtered;
}

Eigen::Vector3d Controller::getZmpFromExternalForces()
{
    bool left_contact = false;
    bool right_contact = false;
    
    Eigen::VectorXd grf_L = mLeftFoot->getConstraintImpulse();
    Eigen::VectorXd grf_R = mRightFoot->getConstraintImpulse();
    Eigen::Vector3d left_cop, right_cop;

    if(abs(grf_L[5]) > 0.01){
        left_cop << -grf_L(1)/grf_L(5), grf_L(0)/grf_L(5), 0.0;
        left_cop = mLeftFoot->getWorldTransform().translation() + mLeftFoot->getWorldTransform().rotation()*left_cop;
        left_contact = true;
    }

    if(abs(grf_R[5]) > 0.01){
        right_cop << -grf_R(1)/grf_R(5), grf_R(0)/grf_R(5), 0.0;
        right_cop = mRightFoot->getWorldTransform().translation() + mRightFoot->getWorldTransform().rotation()*right_cop;
        right_contact = true;
    }

    if(left_contact && right_contact){
        return Eigen::Vector3d((left_cop(0)*grf_L[5] + right_cop(0)*grf_R[5])/(grf_L[5] + grf_R[5]),
                               (left_cop(1)*grf_L[5] + right_cop(1)*grf_R[5])/(grf_L[5] + grf_R[5]),
		               0.0);
    }else if(left_contact){
        return Eigen::Vector3d(left_cop(0), left_cop(1), 0.0);
    }else if(right_contact){
        return Eigen::Vector3d(right_cop(0), right_cop(1), 0.0);
    }else{
        return current.zmpPos;
    }
}

void Controller::setParamsFromFile(const std::string & filepath)  {
  toml::table config = toml::parse_file(filepath);

  vref_const.x = config.at_path("vel[0]").as_floating_point()->get();
  vref_const.y = config.at_path("vel[1]").as_floating_point()->get();
  vref_const.omega = config.at_path("vel[2]").as_floating_point()->get();

  num_steps = config.at_path("num_steps").as_integer()->get();

  push_duration_ = config.at_path("push.duration").as_integer()->get();
  push_timing_ = config.at_path("push.timing").as_integer()->get();
  push_magnitude_ = config.at_path("push.magnitude").as_integer()->get();
  push_angle_ = config.at_path("push.angle").as_integer()->get();

  cov_meas_pos = config.at_path("kalman.cov_meas_pos").as_floating_point()->get();
  cov_meas_vel = config.at_path("kalman.cov_meas_vel").as_floating_point()->get();
  cov_mod_pos = config.at_path("kalman.cov_mod_pos").as_floating_point()->get();
  cov_mod_vel = config.at_path("kalman.cov_mod_vel").as_floating_point()->get();

}

void Controller::setInitialConfiguration() {

  // floating base
  mRobot->setPosition(0, 0.0 );
  mRobot->setPosition(1, -0*M_PI/180.0);
  mRobot->setPosition(2, 0.0 );
  mRobot->setPosition(3, 0.035502257 + 0.0);
  mRobot->setPosition(4, -0.0);
  mRobot->setPosition(5, 0.751 +0.00138 - 0.000005);

  // right leg
  mRobot->setPosition(mRobot->getDof("R_HIP_Y")->getIndexInSkeleton(), 0.0 );
  mRobot->setPosition(mRobot->getDof("R_HIP_R")->getIndexInSkeleton(), -3*M_PI/180 );
  mRobot->setPosition(mRobot->getDof("R_HIP_P")->getIndexInSkeleton(), -25*M_PI/180 );
  mRobot->setPosition(mRobot->getDof("R_KNEE_P")->getIndexInSkeleton(), 50*M_PI/180 );
  mRobot->setPosition(mRobot->getDof("R_ANKLE_P")->getIndexInSkeleton(), -26*M_PI/180+0.0175 );
  mRobot->setPosition(mRobot->getDof("R_ANKLE_R")->getIndexInSkeleton(), 4*M_PI/180-0.01745);

  // left leg
  mRobot->setPosition(mRobot->getDof("L_HIP_Y")->getIndexInSkeleton(), 0.0 );
  mRobot->setPosition(mRobot->getDof("L_HIP_R")->getIndexInSkeleton(), 3*M_PI/180 );
  mRobot->setPosition(mRobot->getDof("L_HIP_P")->getIndexInSkeleton(), -25*M_PI/180 );
  mRobot->setPosition(mRobot->getDof("L_KNEE_P")->getIndexInSkeleton(), 50*M_PI/180 );
  mRobot->setPosition(mRobot->getDof("L_ANKLE_P")->getIndexInSkeleton(), -26*M_PI/180+0.0175 );
  mRobot->setPosition(mRobot->getDof("L_ANKLE_R")->getIndexInSkeleton(), -4*M_PI/180+0.01745);

  // right arm
  mRobot->setPosition(mRobot->getDof("R_SHOULDER_P")->getIndexInSkeleton(), (4)*M_PI/180 );
  mRobot->setPosition(mRobot->getDof("R_SHOULDER_R")->getIndexInSkeleton(), -8*M_PI/180  );
  mRobot->setPosition(mRobot->getDof("R_SHOULDER_Y")->getIndexInSkeleton(), 0 );
  mRobot->setPosition(mRobot->getDof("R_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 );

  // left arm
  mRobot->setPosition(mRobot->getDof("L_SHOULDER_P")->getIndexInSkeleton(), (4)*M_PI/180  );
  mRobot->setPosition(mRobot->getDof("L_SHOULDER_R")->getIndexInSkeleton(), 8*M_PI/180  );
  mRobot->setPosition(mRobot->getDof("L_SHOULDER_Y")->getIndexInSkeleton(), 0 );
  mRobot->setPosition(mRobot->getDof("L_ELBOW_P")->getIndexInSkeleton(), -25*M_PI/180 );

}



void Controller::draw() {

  mRobot->getBodyNode("base_link")->removeAllShapeNodesWith<dart::dynamics::VisualAspect>();

  if (walk_state_.global_iter >= push_timing_ && walk_state_.global_iter < push_timing_ + push_duration_ && push_magnitude_ > 1e-6) {
        dart::dynamics::ArrowShape::Properties arrow_properties;
        arrow_properties.mRadius = 0.05;
        Eigen::Vector3d head_pos = mTorso->getTransform().translation() - mRobot->getBodyNode("base_link")->getTransform().translation();
        Eigen::Vector3d tail_pos = head_pos - 0.005*push_magnitude_* Eigen::Vector3d(-cos(push_angle_/180.0*M_PI),-sin(push_angle_/180.0*M_PI),0.0);
        std::shared_ptr<dart::dynamics::ArrowShape> mArrow = std::shared_ptr<dart::dynamics::ArrowShape>(new dart::dynamics::ArrowShape(tail_pos, head_pos, arrow_properties, dart::Color::Red(1.0)));
        mRobot->getBodyNode("base_link")->createShapeNodeWith<dart::dynamics::VisualAspect>(mArrow, "push_arrow");
  }


/*  int shapeNodesCounterStart = (mGround->getBodyNode("ground_link")->getShapeNodesWith<dart::dynamics::VisualAspect>()).size();
  for (int i = 0; i < footstep_planner->getSize()-10; ++i) {
        auto footstepShapeNode = (mGround->getBodyNode("ground_link")->getShapeNodesWith<dart::dynamics::VisualAspect>()).at(shapeNodesCounterStart + i);
        footstepShapeNode->setRelativeTranslation(footstep_planner->getFootstepPosition(i));
        Eigen::Matrix3d RotZ = rotz(footstepPlan->getFootstepOrientation(i));
        footstepShapeNode->setRelativeRotation(RotZ);
  }*/
}