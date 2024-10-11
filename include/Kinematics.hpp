#pragma once

#include <pinocchio/algorithm/model.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/centroidal.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <dart/dart.hpp>
//#include <dart/utils/utils.hpp>

#include "parameters.hpp"
#include "types.hpp"
#include "utils.hpp"



class Kinematics {

 public:
  Kinematics(const dart::dynamics::SkeletonPtr & robot, const std::string & urdf_filename);

  void updateKinematics(const Eigen::VectorXd & config, const Eigen::VectorXd & joint_vel, Foot predSF);

  Eigen::MatrixXd getTaskJacobian(Task task) const;
  Eigen::MatrixXd getTaskJacobianTimeDeriv(Task task) const;
  Eigen::VectorXd getTaskPosition(Task task) const;
  Eigen::VectorXd getTaskVelocity(Task task) const;

  Eigen::VectorXd predictNextConfig();
  Eigen::VectorXd predictNextConfig(const Eigen::VectorXd & qdot);

  Eigen::MatrixXd getCentroidalJacobian();

  Eigen::VectorXd coords2joints(const Eigen::VectorXd & c);
  Eigen::MatrixXd coords2joints(const Eigen::MatrixXd & A);
  Eigen::VectorXd joints2coords(const Eigen::VectorXd & q);

  static Eigen::VectorXd D2PFloatingBase(const Eigen::Vector6d & dart_pose);

  int getDofs() const;
  int getDofsFullModel() const;


  void debug();

 private:


  int n_joints;
  int n_joints_full;
  std::vector<int> dof_indices; // Converts between the joint list and the dart configuration

  pinocchio::Model robot_model;
  pinocchio::Data current_data;
  Eigen::VectorXd current_config;
  Eigen::VectorXd current_velocity;

  int r_sole_ID;
  int l_sole_ID;
  int torso_ID;
  int CoM_ID;
  int r_hand_ID, l_hand_ID;


  State task_positions_;

  Eigen::MatrixXd J_swing_foot_;
  Eigen::MatrixXd J_CoM_;
  Eigen::MatrixXd J_torso_;
  Eigen::MatrixXd J_supp2fb_;
  Eigen::MatrixXd J_left_foot_;
  Eigen::MatrixXd J_right_foot_;
  Eigen::MatrixXd J_CMM_;
  Eigen::MatrixXd J_right_hand_, J_left_hand_;

  Eigen::MatrixXd Jdot_swing_foot_;
  Eigen::MatrixXd Jdot_CoM_;
  Eigen::MatrixXd Jdot_torso_;
  Eigen::MatrixXd Jdot_supp2fb_;
  Eigen::MatrixXd Jdot_left_foot_;
  Eigen::MatrixXd Jdot_right_foot_;
  Eigen::MatrixXd Jdot_CMM_;
  Eigen::MatrixXd Jdot_right_hand_, Jdot_left_hand_;

  void computeJacobians(Foot support_foot);
  void computeTaskPositions(Foot support_foot);


  pinocchio::Model buildModelFromJoints(const std::string & urdf_filename, const std::vector<std::string> & joint_names);
  Eigen::VectorXd D2PInitialConfig(const pinocchio::Model & model);


};

