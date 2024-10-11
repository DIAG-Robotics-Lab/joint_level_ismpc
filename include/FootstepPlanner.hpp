#pragma once

#include <iostream>
#include <vector>
#include <fstream>

#include <Eigen/Core>
#include "labrob_qpsolvers/qpsolvers.hpp"

#include "types.hpp"
#include "parameters.hpp"
#include "tools/toml.hpp"



class FootstepPlanner {
 public:
  FootstepPlanner();
  ~FootstepPlanner();
  void plan(const std::vector<Vref> & vrefSequence, const Eigen::Vector6d& initialLeftFoot, const Eigen::Vector6d& initialRightFoot, Foot firstSupportFoot, int current_time);
  Eigen::VectorXd getFootstep(int num);
  Eigen::VectorXd getFootstepPosition(int num);
  Eigen::VectorXd getFootstepPose(int num);
  double getFootstepOrientation(int num);
  int getFootstepStartTiming(int num);
  int getFootstepEndTiming(int num);
  int getFootstepDuration(int num);
  int getFootstepIndexAtTime(int time);
  int getMode(int num);
  int getSize();
  std::vector<Eigen::VectorXd> getPlan();

  Eigen::Vector6d getSwingPos(int num);
  Eigen::Vector6d getSwingVel(int num);
  std::vector<Eigen::Vector6d> getSwingPos();
  std::vector<Eigen::Vector6d> getSwingVel();

  Eigen::Vector2d getZMPPos(int num);
  Eigen::Vector2d getZMPVel(int num);
  std::vector<Eigen::Vector2d> getZMPPos();
  std::vector<Eigen::Vector2d> getZMPVel();

  double getTorsoOrientation(int num);
  double getTorsoOrientationVel(int num);
  std::vector<double> getTorsoOrientation();
  std::vector<double> getTorsoOrientationVel();

  bool isSupportFootLeft(int num);

  void computeSwingTrajectory();
  void computeZMPTrajectory();
  void computeModeSequence();

  std::vector<Eigen::VectorXd> footstepPlan;

 private:

  double swing_height;
  double ell;
  double coronalDeviationMax;
  double sagittalDeviationMax;

  std::vector<Eigen::Vector6d> swing_trajectory_pos_;
  std::vector<Eigen::Vector6d> swing_trajectory_vel_;
  std::vector<Eigen::Vector2d> ZMP_trajectory_pos_;
  std::vector<Eigen::Vector2d> ZMP_trajectory_vel_;
  std::vector<double> torso_orientation_;
  std::vector<double> torso_orientation_vel_;
  Eigen::Vector6d initial_left_foot_;
  Eigen::Vector6d initial_right_foot_;
  Foot firstSupportFoot_;
};
