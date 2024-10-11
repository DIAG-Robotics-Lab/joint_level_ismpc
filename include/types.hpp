#pragma once

#include <Eigen/Core>
#include <vector>
#include "utils.hpp"


enum class Task {COM_POSITION, 
                TORSO_POSE, 
                TORSO_ORIENTATION, 
                SWING_FOOT_POSE, 
                SUPPORT_FOOT_POSE, 
                LEFT_FOOT_POSE, 
                RIGHT_FOOT_POSE, 
                CMM, 
                RIGHT_HAND_POSE, 
                LEFT_HAND_POSE,
                HIP_POSE};

enum class Foot {LEFT, RIGHT};

enum class WalkMode {STANDING, WALKING};

inline void switchSupportFoot(Foot& sf) {
  if(sf == Foot::RIGHT) sf = Foot::LEFT;
  else sf = Foot::RIGHT;
}
inline std::string supportFoot2String(const Foot sf) {
  if(sf == Foot::RIGHT) return "RIGHT";
  else return "LEFT";
}

template <typename VectorEigen>
struct Predictions {
  std::vector<Eigen::MatrixXd> jacobians;
  std::vector<Eigen::MatrixXd> jacobians_dot;
  std::vector<VectorEigen> pos;
  inline void clear() {
    jacobians.clear();
    jacobians_dot.clear();
    pos.clear();
  }
};

template <typename VectorEigen>
struct References {
  std::vector<VectorEigen> pos;
  std::vector<VectorEigen> vel;
  inline void clear() {
    pos.clear();
    vel.clear();
  }
};

struct EndEffector {
  Eigen::Vector3d pos, vel, ang_pos, ang_vel;
};

struct State {
  EndEffector CoM, leftFoot, rightFoot, leftHand, rightHand;
  Eigen::Vector3d zmpPos;
  Foot support_foot;

  inline Eigen::Vector6d getCoMPose() const {
    Eigen::Vector6d CoMPose;
    CoMPose << CoM.pos, CoM.ang_pos;
    return CoMPose;
  }

  inline Eigen::Vector6d getCoMVelocity() const {
    Eigen::Vector6d CoMVelocity;
    CoMVelocity <<  CoM.vel, CoM.ang_vel;
    return CoMVelocity;
  }

  inline Eigen::Vector6d getSupportFootPose() const {
    Eigen::Vector6d sfPose;
    if (support_foot == Foot::LEFT) sfPose << leftFoot.pos, leftFoot.ang_pos;
    else sfPose << rightFoot.pos, rightFoot.ang_pos;
    return sfPose;
  }

  inline Eigen::Vector6d getSwingFootPose() const {
    Eigen::Vector6d sfPose;
    if (support_foot == Foot::RIGHT) sfPose << leftFoot.pos, leftFoot.ang_pos;
    else sfPose << rightFoot.pos, rightFoot.ang_pos;
    return sfPose;
  }

  inline Eigen::Vector6d getSwingFootVelocity() const {
    Eigen::Vector6d sfVelocity;
    if (support_foot == Foot::RIGHT) sfVelocity << leftFoot.vel, leftFoot.ang_vel;
    else sfVelocity <<  rightFoot.vel, rightFoot.ang_vel;
    return sfVelocity;
  }


  inline void print() {
    std::cout << "CoM position: " << CoM.pos;
    std::cout << "left foot position: " << leftFoot.pos;
    std::cout << "right foot position: " << rightFoot.pos;
  }
};

struct Data {
  References<Eigen::Vector2d> zmp_refs;
  std::vector<Eigen::Vector2d> zmp_prev;

  std::vector<Eigen::VectorXd> q;
  std::vector<Eigen::VectorXd> qdot;

  References<Eigen::Vector3d> CoM_refs;
  References<Eigen::Vector3d> torso_orientation_refs;
//  References<Eigen::Vector6d> swing_foot_refs;
  References<Eigen::Vector6d> left_foot_refs;
  References<Eigen::Vector6d> right_foot_refs;

  References<Eigen::Vector6d> right_hand_refs, left_hand_refs, hip_refs;

  Predictions<Eigen::Vector3d> centroidal_momentum_preds;
  Predictions<Eigen::Vector3d> CoM_preds;
  Predictions<Eigen::Vector3d> torso_orientation_preds;
//  Predictions<Eigen::Vector6d> swing_foot_preds;
  Predictions<Eigen::Vector6d> left_foot_preds;
  Predictions<Eigen::Vector6d> right_foot_preds;
//  Predictions<Eigen::Vector6d> torso_pose_preds;

  Predictions<Eigen::Vector6d> right_hand_preds, left_hand_preds, hip_preds;

  Eigen::Vector3d filtered_zmp_pos;

  Eigen::VectorXd current_q;
  Eigen::VectorXd q0;
  Eigen::VectorXd current_qdot;

//  Eigen::Vector3d current_com_pos;
//  Eigen::Vector3d current_com_vel;

  State current_state;

  inline void clearReferences() {
    zmp_refs.clear();
    zmp_prev.clear();

    CoM_refs.clear();
    torso_orientation_refs.clear();

    left_foot_refs.clear();
    right_foot_refs.clear();

    right_hand_refs.clear();
    left_hand_refs.clear();
    hip_refs.clear();
  }

  inline void clearAll() {
    zmp_refs.clear();
    zmp_prev.clear();
    q.clear();
    qdot.clear();
    centroidal_momentum_preds.clear();

    CoM_refs.clear();
    torso_orientation_refs.clear();
  //  swing_foot_refs.clear();
    left_foot_refs.clear();
    right_foot_refs.clear();

    CoM_preds.clear();
    torso_orientation_preds.clear();
  //  swing_foot_preds.clear();
    right_foot_preds.clear();
    left_foot_preds.clear();
//    torso_pose_preds.clear();

  //  support_foot_preds.clear();
    right_hand_refs.clear();
    left_hand_refs.clear();
    hip_refs.clear();
    right_hand_preds.clear();
    left_hand_preds.clear();
    hip_preds.clear();
  }
};



struct WalkState {
  Foot support_foot;
  int global_iter, mpc_iter, footstep_idx;
  std::deque<Foot> predicted_support_foot;
  WalkMode mode;
};

struct Gains {
  double CoM_height = 0;
  Eigen::Vector3d torso_orient = Eigen::Vector3d::Zero();
  Eigen::Vector6d swing_pose = Eigen::Vector6d::Zero();
};

struct Vref {
  Vref() {}
  Vref(double _x, double _y, double _omega) : x(_x), y(_y), omega(_omega) {}

  double x = 0;
  double y = 0;
  double omega = 0;
};


