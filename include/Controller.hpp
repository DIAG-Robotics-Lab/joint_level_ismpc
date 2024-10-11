#pragma once

#include <math.h>
#include <fstream>
#include <string>

#include "Kinematics.hpp"
#include "JointLevelISMPC.hpp"
#include "FootstepPlanner.hpp"

class Controller
{
public:
  Controller(dart::dynamics::SkeletonPtr _robot, Kinematics _robot_kinematics);
  virtual ~Controller();

  void update();

  int num_steps;
  // Contains all the timings
  WalkState walk_state_;


private:
  Kinematics robotKinematics;

  void draw();

  Data mpc_data_;

  dart::dynamics::SkeletonPtr mRobot;
  dart::dynamics::BodyNode* mTorso;
  dart::dynamics::BodyNode* mLeftFoot;
  dart::dynamics::BodyNode* mRightFoot;

  Eigen::Vector3d zmp_buffer_box;



  int n_joints;
  int n_joints_full;

  int push_duration_, push_timing_, push_magnitude_, push_angle_;


  std::shared_ptr<FootstepPlanner> footstep_planner;

  std::unique_ptr<JointLevelISMPC> mpc;

  Eigen::MatrixXd QAcc;

  State desired;
  State current;
  State initial;

  Vref vref_const;



  Eigen::VectorXd jointPos, jointVel;
  Eigen::VectorXd swing_des;
  Eigen::MatrixXd swing_pin;
  Eigen::MatrixXd predictedJointPos;
  Eigen::MatrixXd predictedCoMPos;

  Eigen::Vector4d feasibility_bounds;
  
  Eigen::Vector3d filteredZmpPos;

  Eigen::MatrixXd com_jacobian;

  Eigen::Vector6d right_hand_pos_log, right_hand_pos_ref_log;

  utils::Logger logger_;

  Eigen::Vector2d ComAccXY = Eigen::Vector2d::Zero();
  Eigen::Vector2d ComAccXYMeas = Eigen::Vector2d::Zero();
  Eigen::Vector2d ZmpRef;

  void planFootsteps();

  void computePredictions();
  void computeReferences();

  void predictLastSupportFoot();

  void setParamsFromFile(const std::string & filepath);

  // Kalman filter
  State updateKF(State filtered, State current, Eigen::Vector2d comAcc);
  Eigen::MatrixXd cov_x = Eigen::Vector2d(0.27, 2.7).asDiagonal();
  Eigen::MatrixXd cov_y = Eigen::Vector2d(0.27, 2.7).asDiagonal();
  double cov_meas_pos, cov_meas_vel, cov_mod_pos, cov_mod_vel;

  Eigen::Vector3d getZmpFromExternalForces();

  void setInitialConfiguration();
  Eigen::VectorXd getCurrentConfiguration();

  Eigen::VectorXd getCurrentVelocity();

  void setLogger() {

      logger_.setLoggingDirectory("./data");

      logger_.add(current.CoM.pos, "ComPosMeas", true, true);
      logger_.add(current.CoM.vel, "ComVelMeas", true, true);

      logger_.add(mpc_data_.current_state.CoM.pos, "ComPosFilt", true, true);
      logger_.add(mpc_data_.current_state.CoM.vel, "ComVelFilt", true, true);

      logger_.add(swing_des, "SwingFootRef", true, true);

      logger_.add(mpc_data_.current_q, "JointPos", true, true);
      logger_.add(mpc_data_.current_q, "JointPos", true, true);

      logger_.add(ComAccXY, "ComAccXY", true, true);
      logger_.add(ComAccXYMeas, "ComAccXYMeas", true, true);

      logger_.add(current.zmpPos, "ZmpPosMeas", true, true);
      logger_.add(current.leftFoot.pos, "LeftFootPos", true, true);
      logger_.add(current.leftFoot.ang_pos, "LeftFootAngPos", true, true);
      logger_.add(current.rightFoot.pos, "RightFootPos", true, true);
      logger_.add(current.rightFoot.ang_pos, "RightFootAngPos", true, true);

      logger_.add(right_hand_pos_log, "RightHandPos", true, true);
      logger_.add(right_hand_pos_ref_log, "RightHandPosRef", true, true);

      logger_.add(ZmpRef, "ZmpRef", true, true);

      logger_.add(cov_x, "CovXKF", false, true);
      logger_.add(cov_y, "CovYKF", false, true);
//
//    logger_.add(new utils::Container<Eigen::VectorXd>(&jointPos, "joint_pos_dart"));
//    logger_.add(new utils::Container<Eigen::VectorXd>(&jointVel, "joint_vel_dart"));
//    logger_.add(new utils::Container<Eigen::MatrixXd>(&predictedJointPos, "joint_pos_pred"));
//
//    logger_.add(new utils::Container<Eigen::VectorXd>(&swing_des, "swing_ref"));
//    logger_.add(new utils::Container<Eigen::MatrixXd>(&swing_pin, "swing_pred"));
//
//    logger_.add(new utils::Container<Eigen::MatrixXd>(&predictedCoMPos, "com_pred_pinocchio"));
//    logger_.add(new utils::Container<Eigen::Vector3d>(&filteredZmpPos, "zmp_filtered"));
//
//    logger_.add(new utils::Container<Eigen::MatrixXd>(&com_jacobian, "com_jacobiankp1"));
//    logger_.add(new utils::Container<Eigen::Vector4d>(&feasibility_bounds, "feasibility_bounds"));
//
//    logger_.add(new utils::Container<Eigen::Vector6d>(&right_hand_pos_log, "right_hand_pos"));
//    logger_.add(new utils::Container<Eigen::Vector6d>(&right_hand_pos_ref_log, "right_hand_pos_ref"));


  };
};
