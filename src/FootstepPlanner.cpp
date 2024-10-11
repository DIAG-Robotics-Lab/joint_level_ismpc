#include "FootstepPlanner.hpp"

FootstepPlanner::FootstepPlanner() {

  toml::table config = toml::parse_file("./simulation_params.toml");

  swing_height = config.at_path("swing_height").as_floating_point()->get();
  ell = config.at_path("ell").as_floating_point()->get();
  sagittalDeviationMax = config.at_path("sagittalDeviationMax").as_floating_point()->get();
  coronalDeviationMax = config.at_path("coronalDeviationMax").as_floating_point()->get();

}

FootstepPlanner::~FootstepPlanner() = default;

void FootstepPlanner::plan(const std::vector<Vref>& vrefSequence,
                           const Eigen::Vector6d& initialLeftFoot,
                           const Eigen::Vector6d& initialRightFoot,
                           Foot firstSupportFoot, int current_time) {

  firstSupportFoot_ = firstSupportFoot;
  // For consistency, we have to switch support foot due to the dummy initial step
  //switchSupportFoot(firstSupportFoot_);
  initial_left_foot_ = initialLeftFoot;
  initial_right_foot_ = initialRightFoot;

  int nSamples = vrefSequence.size()+1;
  std::vector<double> stepStartingSequence; // when each step starts (in MPC samples)

  stepStartingSequence.push_back(0);

  // Integrate reference velocities using omnidirectional model

  Eigen::VectorXd integrated_x(nSamples);
  Eigen::VectorXd integrated_y(nSamples);
  Eigen::VectorXd integrated_theta(nSamples);
  integrated_x(0) = (initialLeftFoot(0) + initialRightFoot(0)) / 2.0;
  integrated_y(0) = (initialLeftFoot(1) + initialRightFoot(1)) / 2.0;
  integrated_theta(0) = (initialLeftFoot(5) + initialRightFoot(5)) / 2.0;

  for (int i = 1; i < nSamples; i++) {
    integrated_theta(i) = integrated_theta(i-1) + vrefSequence.at(i-1).omega * worldTimeStep;
    Eigen::Matrix2d rotationMatrix;
    rotationMatrix << cos(integrated_theta(i)), -sin(integrated_theta(i)), sin(integrated_theta(i)), cos(integrated_theta(i));
    Eigen::Vector2d temp = Eigen::Vector2d(integrated_x(i-1), integrated_y(i-1)) + rotationMatrix * Eigen::Vector2d(vrefSequence.at(i-1).x, vrefSequence.at(i-1).y) * worldTimeStep;
    integrated_x(i) = temp(0);
    integrated_y(i) = temp(1);
  }


  int timingIndex = 0;
  while (timingIndex < nSamples-1) {
    //double velocityMagnitude = sqrt(vrefSequence.at(timingIndex).x * vrefSequence.at(timingIndex).x + vrefSequence.at(timingIndex).y * vrefSequence.at(timingIndex).y);
    double stepDuration = singleSupportDuration + doubleSupportDuration;
//         if (stepDuration >= maxStepDuration) stepDuration = maxStepDuration;
//         if (stepDuration <= minStepDuration) stepDuration = minStepDuration;

    timingIndex = timingIndex + (int)(stepDuration / worldTimeStep);
    stepStartingSequence.push_back(timingIndex);
  }

  int nFootsteps = stepStartingSequence.size();
  Eigen::VectorXd xFootsteps(nFootsteps);
  Eigen::VectorXd yFootsteps(nFootsteps);
  Eigen::VectorXd thetaFootsteps(nFootsteps);

  // Select integrated footstep positions and orientations at planned timings
  Eigen::VectorXd xFootstepsSelected(nFootsteps);
  Eigen::VectorXd yFootstepsSelected(nFootsteps);
  Eigen::VectorXd thetaFootstepsSelected(nFootsteps);

  if (firstSupportFoot == Foot::LEFT) {
    xFootstepsSelected(0) = initialLeftFoot(0);
    yFootstepsSelected(0) = initialLeftFoot(1);
    thetaFootstepsSelected(0) = initialLeftFoot(5);
  } else {
    xFootstepsSelected(0) = initialRightFoot(0);
    yFootstepsSelected(0) = initialRightFoot(1);
    thetaFootstepsSelected(0) = initialRightFoot(5);
  }

  for (int i = 1; i < nFootsteps; i++) {
    thetaFootstepsSelected(i) = integrated_theta((int)stepStartingSequence.at(i));
    Eigen::Matrix2d ortoRotationMatrix;

    double ortoTheta;
    if (firstSupportFoot == Foot::LEFT) ortoTheta = thetaFootstepsSelected(i) + isEven(i) * M_PI/2.0 - isOdd(i) * M_PI/2.0;
    else ortoTheta = thetaFootstepsSelected(i) + isOdd(i) * M_PI/2.0 - isEven(i) * M_PI/2.0;

    ortoRotationMatrix << cos(ortoTheta), -sin(ortoTheta), sin(ortoTheta), cos(ortoTheta);
    Eigen::Vector2d temp = Eigen::Vector2d(integrated_x((int)stepStartingSequence.at(i)), integrated_y((int)stepStartingSequence.at(i))) + ortoRotationMatrix * Eigen::Vector2d(ell/2.0, 0.0);
    xFootstepsSelected(i) = temp(0);
    yFootstepsSelected(i) = temp(1);
    thetaFootstepsSelected(i) = integrated_theta((int)stepStartingSequence.at(i));
  }

  // Plan orientations using QP

  Eigen::MatrixXd differenceMatrix = Eigen::MatrixXd::Identity(nFootsteps, nFootsteps);
  for (int i = 0; i < nFootsteps - 1; i++) {
    differenceMatrix(i+1, i) = - 1;
  }

  Eigen::MatrixXd orientationsH = Eigen::MatrixXd::Identity(nFootsteps, nFootsteps);
  Eigen::VectorXd orientationsF = - thetaFootstepsSelected;
  Eigen::MatrixXd orientationsA = differenceMatrix;
  Eigen::VectorXd orientationsBmax = Eigen::VectorXd::Ones(nFootsteps) * thetaMax;
  Eigen::VectorXd orientationsBmin = - Eigen::VectorXd::Ones(nFootsteps) * thetaMax;

  int numVariables_ = orientationsF.size();
  int numInequalityConstraints_ = orientationsBmax.size();
  auto orientation_QP = std::make_shared<labrob::qpsolvers::QPSolverEigenWrapper<double>>(
      std::make_shared<labrob::qpsolvers::HPIPMQPSolver>(
          numVariables_, 1, numInequalityConstraints_));
  // for now dummy equality constraint since the qp interface requires it
  Eigen::MatrixXd dummy_constraint = 0*Eigen::MatrixXd::Zero(1, numVariables_);
  Eigen::VectorXd dummy_scalar= 0*Eigen::VectorXd::Zero(1);

  orientation_QP->solve(
      orientationsH,
      orientationsF,
      dummy_constraint,
      dummy_scalar,
      orientationsA,
      orientationsBmin,
      orientationsBmax
  );
  thetaFootsteps = orientation_QP->get_solution();//solveQP(orientationsH, orientationsF, orientationsA, orientationsBmin, orientationsBmax);

  for (int i = 0; i < nFootsteps; i++) {
    thetaFootsteps(i) = wrapToPi(thetaFootsteps(i));
  }

  // Plan positions using QP

  Eigen::MatrixXd doubleDifferenceMatrix = Eigen::MatrixXd::Identity(2 * nFootsteps, 2 * nFootsteps);
  for (int i = 0; i < 2 * (nFootsteps - 1); i++) {
    doubleDifferenceMatrix(i+2, i) = - 1;
  }

  Eigen::MatrixXd rotationMatrices = Eigen::MatrixXd::Zero(2 * nFootsteps, 2 * nFootsteps);
  for (int i = 0; i < nFootsteps; i++) {
    Eigen::Matrix2d rotationBlock;
    if (i == 0) {
      rotationBlock = Eigen::Matrix2d::Identity();
    } else {
      rotationBlock << cos(thetaFootsteps(i-1)), -sin(thetaFootsteps(i-1)), sin(thetaFootsteps(i-1)), cos(thetaFootsteps(i-1));
    }

    rotationMatrices.block(2*i, 2*i, 2, 2) = rotationBlock.transpose();
  }

  Eigen::VectorXd ellVector = Eigen::VectorXd::Zero(2 * nFootsteps);
  Eigen::VectorXd integratedFootstepsVector = Eigen::VectorXd::Zero(2 * nFootsteps);
  Eigen::VectorXd boxSizeVector = Eigen::VectorXd::Zero(2 * nFootsteps);
  for (int i = 0; i < nFootsteps; i++) {
    if (firstSupportFoot == Foot::LEFT) ellVector(2*i + 1) = isEven(i) ? ell : -ell;
    else ellVector(2*i + 1) = isOdd(i) ? ell : -ell;
    integratedFootstepsVector(2*i) = xFootstepsSelected(i);
    integratedFootstepsVector(2*i + 1) = yFootstepsSelected(i);
    boxSizeVector(2*i) = sagittalDeviationMax;
    boxSizeVector(2*i + 1) = coronalDeviationMax;
  }

  Eigen::VectorXd initialFootstepVector = Eigen::VectorXd::Zero(2 * nFootsteps);
  initialFootstepVector(0) = integratedFootstepsVector(0);
  initialFootstepVector(1) = -integratedFootstepsVector(1);


  Eigen::MatrixXd positionsH = Eigen::MatrixXd::Identity(2 * nFootsteps, 2 * nFootsteps);
  Eigen::VectorXd positionsF = - integratedFootstepsVector;
  Eigen::MatrixXd positionsA = rotationMatrices * doubleDifferenceMatrix;
  Eigen::VectorXd positionsBmax = ellVector + boxSizeVector + initialFootstepVector;
  Eigen::VectorXd positionsBmin = ellVector - boxSizeVector + initialFootstepVector;

  numVariables_ = positionsF.size();
  numInequalityConstraints_ = positionsBmax.size();
  auto positions_QP = std::make_shared<labrob::qpsolvers::QPSolverEigenWrapper<double>>(
      std::make_shared<labrob::qpsolvers::HPIPMQPSolver>(
          numVariables_, 1, numInequalityConstraints_));

  dummy_constraint = Eigen::MatrixXd::Zero(1, numVariables_);
  dummy_scalar= Eigen::VectorXd::Zero(1);

  positions_QP->solve(
      positionsH,
      positionsF,
      dummy_constraint,
      dummy_scalar,
      positionsA,
      positionsBmin,
      positionsBmax
  );

  Eigen::VectorXd positionsDecisionVariables(2 * nFootsteps);
  positionsDecisionVariables = positions_QP->get_solution();


  for (int i = 0; i < nFootsteps; i++) {
    xFootsteps(i) = integratedFootstepsVector(2*i);
    yFootsteps(i) = integratedFootstepsVector(2*i + 1);
  }

  /*
  // Mode sequence
  std::vector<int> modeSequence;
  modeSequence.push_back(0);
  modeSequence.push_back(0);
  for (int i = 2; i < nFootsteps-1; i++) {
    // if the next footstep is unchanged, the mode is 0
    if (abs(getFootstepPosition(i) - getFootstepPosition(i-2)).norm() < 1e-6) {
      modeSequence.push_back(0);
    } else {
      modeSequence.push_back(1);
    }
  }
  modeSequence.push_back(0);
  */

  // Fill footstep plan

  const double z_offset = initial_right_foot_(2);

  for (int i = 0; i < nFootsteps; i++) {
    Eigen::VectorXd tempFootstep(8);
    tempFootstep << xFootsteps(i), yFootsteps(i), z_offset, 0.0, 0.0, thetaFootstepsSelected(i), stepStartingSequence.at(i), 0;
    footstepPlan.push_back(tempFootstep);
  }
  
  // These are added manually to not make the simulation crash at the end
  footstepPlan.push_back(footstepPlan[footstepPlan.size()-2]);
  footstepPlan.push_back(footstepPlan[footstepPlan.size()-2]);
  footstepPlan.push_back(footstepPlan[footstepPlan.size()-2]);
  footstepPlan.push_back(footstepPlan[footstepPlan.size()-2]);
  footstepPlan.push_back(footstepPlan[footstepPlan.size()-2]);
  footstepPlan.push_back(footstepPlan[footstepPlan.size()-2]);

  computeModeSequence();

  // Hack: add current time to all timings
  /*for (int i = 0; i < footstepPlan.size(); i++) {
    footstepPlan.at(i)(6) += current_time;
  }*/


  // Write to file and plot

  /*
  std::ofstream foutReferenceVelocities(realpath("../data/referenceVelocities.txt", NULL), std::ofstream::trunc);
  for (int i = 0; i < vrefSequence.size(); i++) {
      foutReferenceVelocities << vrefSequence.at(i).x << " " << vrefSequence.at(i).y << " " << vrefSequence.at(i).omega << std::endl;
  }
  foutReferenceVelocities.close();

  std::ofstream foutIntegratedTrajectory(realpath("../data/integratedTrajectory.txt", NULL), std::ofstream::trunc);
  for (int i = 0; i < vrefSequence.size(); i++) {
      foutIntegratedTrajectory << integrated_x(i) << " " << integrated_y(i) << std::endl;
  }
  foutIntegratedTrajectory.close();

  std::ofstream foutSelectedFootsteps(realpath("../data/selectedFootsteps.txt", NULL), std::ofstream::trunc);
  for (int i = 0; i < xFootstepsSelected.size(); i++) {
      foutSelectedFootsteps << xFootstepsSelected(i) << " " << yFootstepsSelected(i) << " " << thetaFootstepsSelected(i) << std::endl;
  }
  foutSelectedFootsteps.close();

  std::ofstream foutFootstepPlan(realpath("../data/footstepPlan.txt", NULL), std::ofstream::trunc);
  for (int i = 0; i < footstepPlan.size(); i++) {
      foutFootstepPlan << footstepPlan.at(i).transpose() << std::endl;
  }
  foutFootstepPlan.close();

  system("gnuplot ../plotters/plotReferenceVelocities");
  */
}

void FootstepPlanner::computeModeSequence() {
  // Mode sequence
  std::vector<int> modeSequence;
  modeSequence.push_back(0);
  modeSequence.push_back(0);
  for (int i = 2; i < footstepPlan.size()-1; i++) {
    // if the next footstep is unchanged, the mode is 0
    if ((getFootstepPosition(i) - getFootstepPosition(i-2)).norm() < 1e-6) {
      modeSequence.push_back(0);
    } else {
      modeSequence.push_back(1);
    }
    
  }
  modeSequence.push_back(0);
  for (int i = 0; i < footstepPlan.size(); i++) footstepPlan.at(i)(7) = modeSequence.at(i);
}

void FootstepPlanner::computeSwingTrajectory() {
  swing_trajectory_pos_.clear();
  swing_trajectory_vel_.clear();
  torso_orientation_.clear();
  torso_orientation_vel_.clear();
  Eigen::VectorXd swing_start;
  if(firstSupportFoot_ == Foot::LEFT) swing_start = initial_right_foot_;
  else swing_start = initial_left_foot_;

  const double z_offset = swing_start(2);

  for(int step = 1; step <getSize(); ++step) {
    // If it's the first step, or we are in double support, keep the foot on the ground
    Eigen::VectorXd swing_end(this->getFootstep(step).head(6));
//    swing_end(2) += z_offset;

    for (int k = 0; k < singleSupportSamples; ++k) {
      double swingFootHeight = z_offset;
      double swingFootVerticalVelocity = 0.0;
      if(getMode(step) == 1) {
        swingFootHeight = z_offset - (4 * swing_height / pow((double) singleSupportSamples, 2)) * k * (k - singleSupportSamples);
        swingFootVerticalVelocity = - (4 * swing_height / pow((double) singleSupportSamples, 2)) * (2 * k - singleSupportSamples)/worldTimeStep;
      }
      Eigen::Vector6d swing_pos_k = Eigen::Vector6d::Zero();
      Eigen::Vector6d swing_vel_k = Eigen::Vector6d::Zero();
      double tau = (k/(double)singleSupportSamples);

      swing_pos_k(0) = swing_start(0) + tau*(swing_end(0)-swing_start(0));
      swing_pos_k(1) = swing_start(1) + tau*(swing_end(1)-swing_start(1));
      swing_pos_k(2) = swingFootHeight;
      swing_pos_k(5) = swing_start(5) + tau*wrapToPi(swing_end(5)-swing_start(5));
      swing_trajectory_pos_.emplace_back(swing_pos_k);

      torso_orientation_.push_back(wrapToPi((this->getFootstep(step-1)(5) + swing_pos_k(5)) / 2.0));

      swing_vel_k(0) = (swing_end(0)-swing_start(0))/((double)singleSupportSamples*worldTimeStep);
      swing_vel_k(1) = (swing_end(1)-swing_start(1))/((double)singleSupportSamples*worldTimeStep);
      swing_vel_k(2) = swingFootVerticalVelocity;
      swing_vel_k(5) = (swing_end(5)-swing_start(5))/((double)singleSupportSamples*worldTimeStep);
      swing_trajectory_vel_.emplace_back(swing_vel_k);

      torso_orientation_vel_.push_back(swing_vel_k(5)/2.0);
    }
    for (int kk = 0; kk < doubleSupportSamples; ++kk) {
      swing_trajectory_pos_.emplace_back(swing_end.head(6));
      swing_trajectory_vel_.emplace_back(Eigen::Vector6d::Zero());

      torso_orientation_.push_back(wrapToPi((this->getFootstep(step-1)(5) + swing_end(5)) / 2.0));
      torso_orientation_vel_.push_back(0.0);
    }

    swing_start = this->getFootstep(step-1);
  }
}

void FootstepPlanner::computeZMPTrajectory() {
  ZMP_trajectory_pos_.clear();
  ZMP_trajectory_vel_.clear();
  Eigen::Vector2d zmp_start = 0.5*(initial_right_foot_.head(2) + initial_left_foot_.head(2));
  for(int step = 1; step<getSize()-1; ++step) {
    Eigen::Vector2d zmp_end;
    if (getMode(step) == 1 || getMode(step+1) == 1) {
      zmp_end = this->getFootstepPosition(step).head(2);
    } else {
      zmp_end = 0.5*(this->getFootstepPosition(step).head(2) + this->getFootstepPosition(step-1).head(2));
    }
    

    for (int k = 0; k < singleSupportSamples; ++k) {
      ZMP_trajectory_pos_.emplace_back(zmp_start);
      ZMP_trajectory_vel_.emplace_back(Eigen::Vector2d::Zero());
    }
    for (int k = 0; k < doubleSupportSamples; ++k) {
      double tau = (k/(double)doubleSupportSamples);
      Eigen::Vector2d zmp_pos_k = Eigen::Vector2d::Zero();
      Eigen::Vector2d zmp_vel_k = Eigen::Vector2d::Zero();

      zmp_pos_k = zmp_start + tau*(zmp_end-zmp_start);
      zmp_vel_k = (zmp_end-zmp_start)/((double)singleSupportSamples*worldTimeStep);
      ZMP_trajectory_pos_.emplace_back(zmp_pos_k);
      ZMP_trajectory_vel_.emplace_back(zmp_vel_k);
    }
  zmp_start = zmp_end;
  }
}

int FootstepPlanner::getSize() { return footstepPlan.size(); }
std::vector<Eigen::VectorXd> FootstepPlanner::getPlan() { return footstepPlan; }

std::vector<Eigen::Vector6d> FootstepPlanner::getSwingPos() { return swing_trajectory_pos_; }
std::vector<Eigen::Vector6d> FootstepPlanner::getSwingVel() { return swing_trajectory_vel_; }
Eigen::Vector6d FootstepPlanner::getSwingPos(int num) {
  try {
    return swing_trajectory_pos_.at(num);
  }
  catch (const std::out_of_range& e) {
    return swing_trajectory_pos_.back();
  }
}
Eigen::Vector6d FootstepPlanner::getSwingVel(int num) {
  return swing_trajectory_vel_.at(num);
}
/*Eigen::Vector6d FootstepPlanner::getSwingVel(int num) {
  try {
    return swing_trajectory_vel_.at(num);
  }
  catch (const std::out_of_range& e) {
    return Eigen::Vector6d::Zero();
  }
}*/

std::vector<double> FootstepPlanner::getTorsoOrientation() { return torso_orientation_; }
std::vector<double> FootstepPlanner::getTorsoOrientationVel() { return torso_orientation_vel_; }
double FootstepPlanner::getTorsoOrientation(int num) {
  try {
    return torso_orientation_.at(num);
  }
  catch (const std::out_of_range& e) {
    return torso_orientation_.back();
  }
}
double FootstepPlanner::getTorsoOrientationVel(int num) {
  try {
    return torso_orientation_vel_.at(num);
  }
  catch (const std::out_of_range& e) {
    return 0.0;
  }
}

std::vector<Eigen::Vector2d> FootstepPlanner::getZMPPos() { return ZMP_trajectory_pos_; }
std::vector<Eigen::Vector2d> FootstepPlanner::getZMPVel() { return ZMP_trajectory_vel_; }
Eigen::Vector2d FootstepPlanner::getZMPPos(int num) {
  try{
    return ZMP_trajectory_pos_.at(num);
  }
  catch (const std::out_of_range& e) {
    return ZMP_trajectory_pos_.back();
  }
}
Eigen::Vector2d FootstepPlanner::getZMPVel(int num) {
  try {
    return ZMP_trajectory_vel_.at(num);
  }
  catch (const std::out_of_range& e) {
    return Eigen::Vector2d::Zero();
  }
}

Eigen::VectorXd FootstepPlanner::getFootstep(int num) { return footstepPlan.at(num); }
Eigen::VectorXd FootstepPlanner::getFootstepPosition(int num) { return footstepPlan.at(num).head(3); }
Eigen::VectorXd FootstepPlanner::getFootstepPose(int num) { return footstepPlan.at(num).head(6); }
double FootstepPlanner::getFootstepOrientation(int num) { return footstepPlan.at(num)(5); }
int FootstepPlanner::getFootstepStartTiming(int num) { return (int)footstepPlan.at(num)(6); }
int FootstepPlanner::getFootstepEndTiming(int num) { return (num <= footstepPlan.size()) ? (int)footstepPlan.at(num + 1)(6) : -1; }
int FootstepPlanner::getFootstepDuration(int num) { return (num <= footstepPlan.size()) ? (int)(footstepPlan.at(num + 1)(6)-footstepPlan.at(num)(6)) : -1; }
int FootstepPlanner::getFootstepIndexAtTime(int time) {
  int footstepIndex = 0;
  while (getFootstepEndTiming(footstepIndex) < time and footstepIndex < footstepPlan.size()) footstepIndex++;
  return footstepIndex;
}

int FootstepPlanner::getMode(int num) { return (int)footstepPlan.at(num)(7); }

bool FootstepPlanner::isSupportFootLeft(int num) {
  if (firstSupportFoot_ == Foot::LEFT) {
    if (num % 2 == 0) return true;
    else return false;
  } else {
    if (num % 2 == 0) return false;
    else return true;
  }
}
