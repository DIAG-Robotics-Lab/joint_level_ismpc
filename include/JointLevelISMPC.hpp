#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <deque>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <labrob_qpsolvers/qpsolvers.hpp>

#include "parameters.hpp"
#include "types.hpp"
#include "utils.hpp"
#include "tools/toml.hpp"

#include "cpp-utils/simulation.h"

#include "FootstepPlanner.hpp"



class JointLevelISMPC {

 public:
  JointLevelISMPC(int n_dof, const Eigen::MatrixXd & joint_lims, const Eigen::VectorXd & q_neutral_tot);
  void prepareMPC(const Data & data, const WalkState & walkstate);
  Eigen::MatrixXd computeQAcc(const Data & data, const WalkState & walkstate);
  std::tuple<double, double, double, double> getFeasibilityRegionBounds(const Data & data);

 private:

  int simulation_outcome = 0;
  int n_dof_;
  int nq_;
  int n_joints_;
  int n_zmp_;
  int numVariables_;
  int numEqualityConstraints_;
  int numInequalityConstraints_;

  int n_joint_pos_constr_;

  Gains reg_gains;

  // Weights
  Eigen::Vector3d avrg_ang_vel_weight;
  Eigen::Vector3d Ldot_weight;
  double joint_pos_weight;
  double torso_weight;
  double swing_weight;
  double zmp_weight;
  double feas_weight;
  double CoM_height_weight;
  double qdot_weight;
  double qddot_weight;
  double footstep_weight;
  double double_support_weight;
  double joint_vel_weight;

  std::deque<Eigen::Vector2d> previous_vel_deque_;


  // step_iter_ tells at which point of the step we are
  //int step_iter_ = 0;

  //Eigen::MatrixXd support_transition_matrix_;

  // Neutral configuration of the robot
  Eigen::VectorXd q_neutral_;

  Eigen::MatrixXd rk_;
  Eigen::MatrixXd rk_prev_; // used for the truncated tail
  Eigen::MatrixXd summation_;
  Eigen::MatrixXd summation_pred_;

  // Stability constraint
  Eigen::MatrixXd A_stab_;
  Eigen::VectorXd b_stab_;

  // ZMP constraint
  Eigen::MatrixXd A_zmp_;
  Eigen::VectorXd b_zmp_L_;
  Eigen::VectorXd b_zmp_R_;
  Eigen::VectorXd mk;

  // Joint position constraint
  Eigen::MatrixXd A_joint_pos_;
  Eigen::VectorXd b_joint_pos_L_;
  Eigen::VectorXd b_joint_pos_R_;


  Eigen::VectorXd kD_lim_L_;
  Eigen::VectorXd kD_lim_R_;


  Eigen::VectorXd current_vel_vector_;
  Eigen::Vector2d angmom_in_zmp_;

  Eigen::MatrixXd P_task_gen_;
  Eigen::MatrixXd p_task_gen_;
  Eigen::MatrixXd V_task_gen_;
  Eigen::MatrixXd v_task_gen_;

  // Final QP matrices
  Eigen::MatrixXd costFunctionH;
  Eigen::VectorXd costFunctionF;
  Eigen::MatrixXd A_eq;
  Eigen::VectorXd b_eq;
  Eigen::MatrixXd A_ineq;
  Eigen::VectorXd b_ineq_L;
  Eigen::VectorXd b_ineq_R;

  Eigen::MatrixXd P_xyCoM, V_xyCoM, P_left_foot, P_right_foot;
  Eigen::VectorXd p_xyCoM, v_xyCoM, p_left_foot, p_right_foot;

  //Logging
  utils::Logger logger_;
  Eigen::VectorXd mk_log_, p_xyCoM_log_;

  Eigen::MatrixXd Mk_log_, P_xyCoM_log_;

  Eigen::VectorXd pred_left_pos, pred_right_pos,  pred_zmp;

  Eigen::VectorXd com_pred_jacobians;

  Eigen::VectorXd xu_star;
  Eigen::VectorXd com_dot_pred_jacobians;

  Eigen::MatrixXd alpha_matrix_;

  std::vector<double> alpha_;

  Eigen::MatrixXd selectX_, selectY_;

  Eigen::Vector2d COMStateToFeasibilityRegion_;

  std::vector<int> joint_regulation_weights_;




  //
  Eigen::MatrixXd joint_pos_matrix;
  Eigen::VectorXd joint_pos_vector;

  Eigen::MatrixXd CoM_height_gains;
  Eigen::MatrixXd torso_orient_gains;
  Eigen::MatrixXd swing_pose_gains;

  // Kinematics
  Eigen::MatrixXd J_CoM;
  Eigen::VectorXd z_CoM_pred_pos;
  Eigen::VectorXd z_CoM_ref_pos;
  Eigen::VectorXd z_CoM_ref_vel;
  Eigen::VectorXd avrg_ang_velocity;


  std::shared_ptr<labrob::qpsolvers::QPSolverEigenWrapper<double>> qp_solver_ptr_;


  Eigen::MatrixXd joints2zmp(const std::vector<Eigen::MatrixXd> & predJs, const Eigen::MatrixXd & P_xyCoM) const;

  Eigen::MatrixXd joints2zmp_angular_momentum(const std::vector<Eigen::MatrixXd> & CoM_jacobians,
                                              const std::vector<Eigen::MatrixXd> & CMM,
                                              const Eigen::MatrixXd & P_xyCoM) const;


  void buildTaskPosMatrix(const std::vector<Eigen::MatrixXd> & jacobians,
                          const std::vector<Eigen::MatrixXd> & jacobians_dot,
                          const std::vector<Eigen::VectorXd> & qdot,
                          const Eigen::VectorXd & task_pos,
                          const Eigen::VectorXd & task_vel,
                          Eigen::MatrixXd & P_task,
                          Eigen::VectorXd & p_task);

  void buildTaskVelMatrix(const std::vector<Eigen::MatrixXd> & jacobians,
                          const std::vector<Eigen::MatrixXd> & jacobians_dot,
                          const std::vector<Eigen::VectorXd> & qdot,
                          const Eigen::VectorXd & task_vel,
                          Eigen::MatrixXd & V_task,
                          Eigen::VectorXd & v_task);

  void computeStabilityConstraint(const Eigen::VectorXd & pk,
                                  const Eigen::MatrixXd & PN,
                                  const Eigen::VectorXd & Xz_ant);

  void computeZMPConstraint(const Eigen::VectorXd & mk,
                            const Eigen::MatrixXd & Mk,
                            const std::vector<Eigen::Vector2d> & zmp_pos,
                            const Data & data);

  void computeAutoZMPConstraint(const Eigen::VectorXd & xck,
                                const Eigen::VectorXd & xrk,
                                const Eigen::VectorXd & xlk,
                                const Eigen::MatrixXd & Mk,
                                const Eigen::MatrixXd & Mr,
                                const Eigen::MatrixXd & Ml,
                                const Data & data);
  
  Eigen::Vector2d COMStateToFeasibilityRegion() {
    Eigen::Vector2d ref;
    ref << (pow(timeStep,20)*pow(omega,20))/1024 + (pow(timeStep,19)*pow(omega,19))/512 + (7*pow(timeStep,18)*pow(omega,18))/128 + (13*pow(timeStep,17)*pow(omega,17))/128 + (317*pow(timeStep,16)*pow(omega,16))/256 + (269*pow(timeStep,15)*pow(omega,15))/128 + (935*pow(timeStep,14)*pow(omega,14))/64 + (355*pow(timeStep,13)*pow(omega,13))/16 + (3095*pow(timeStep,12)*pow(omega,12))/32 + (2045*pow(timeStep,11)*pow(omega,11))/16 + 362*pow(timeStep,10)*pow(omega,10) + (1603*pow(timeStep,9)*pow(omega,9))/4 + (5913*pow(timeStep,8)*pow(omega,8))/8 + (2601*pow(timeStep,7)*pow(omega,7))/4 + (1521*pow(timeStep,6)*pow(omega,6))/2 + 492*pow(timeStep,5)*pow(omega,5) + (1365*pow(timeStep,4)*pow(omega,4))/4 + (285*pow(timeStep,3)*pow(omega,3))/2 + 50*pow(timeStep,2)*pow(omega,2) + 10*timeStep*omega + 1,
        (pow(timeStep,19)*pow(omega,19) + 2*pow(timeStep,18)*pow(omega,18) + 52*pow(timeStep,17)*pow(omega,17) + 96*pow(timeStep,16)*pow(omega,16) + 1076*pow(timeStep,15)*pow(omega,15) + 1800*pow(timeStep,14)*pow(omega,14) + 11360*pow(timeStep,13)*pow(omega,13) + 16800*pow(timeStep,12)*pow(omega,12) + 65440*pow(timeStep,11)*pow(omega,11) + 82752*pow(timeStep,10)*pow(omega,10) + 205184*pow(timeStep,9)*pow(omega,9) + 211968*pow(timeStep,8)*pow(omega,8) + 332928*pow(timeStep,7)*pow(omega,7) + 263424*pow(timeStep,6)*pow(omega,6) + 251904*pow(timeStep,5)*pow(omega,5) + 138240*pow(timeStep,4)*pow(omega,4) + 72960*pow(timeStep,3)*pow(omega,3) + 23040*pow(timeStep,2)*pow(omega,2) + 5120*timeStep*omega + 512)/(512*omega);
    return ref;
  }


  Eigen::Vector2d DesFeasibilityFromXzBounds(const Eigen::VectorXd & zmp_pos_ref) {

    Eigen::VectorXd XYz_max = zmp_pos_ref + kD_lim_R_;
    Eigen::VectorXd XYz_min = zmp_pos_ref + kD_lim_L_;

    Eigen::VectorXd Xz_max = selectX_*XYz_max;
    Eigen::VectorXd Xz_min = selectX_*XYz_min;
    Eigen::VectorXd Yz_max = selectY_*XYz_max;
    Eigen::VectorXd Yz_min = selectY_*XYz_min;

    Eigen::VectorXd Pc_times_invM(10);
    Pc_times_invM << -(timeStep*omega*(pow(timeStep,19)*pow(omega,19) + 2*pow(timeStep,18)*pow(omega,18) + 54*pow(timeStep,17)*pow(omega,17) + 100*pow(timeStep,16)*pow(omega,16) + 1168*pow(timeStep,15)*pow(omega,15) + 1968*pow(timeStep,14)*pow(omega,14) + 12992*pow(timeStep,13)*pow(omega,13) + 19456*pow(timeStep,12)*pow(omega,12) + 79584*pow(timeStep,11)*pow(omega,11) + 102592*pow(timeStep,10)*pow(omega,10) + 268096*pow(timeStep,9)*pow(omega,9) + 284544*pow(timeStep,8)*pow(omega,8) + 472320*pow(timeStep,7)*pow(omega,7) + 387072*pow(timeStep,6)*pow(omega,6) + 391680*pow(timeStep,5)*pow(omega,5) + 224256*pow(timeStep,4)*pow(omega,4) + 125184*pow(timeStep,3)*pow(omega,3) + 41472*pow(timeStep,2)*pow(omega,2) + 9728*timeStep*omega + 1024))/1024,
        -(timeStep*omega*(pow(timeStep,17)*pow(omega,17) + 2*pow(timeStep,16)*pow(omega,16) + 48*pow(timeStep,15)*pow(omega,15) + 88*pow(timeStep,14)*pow(omega,14) + 896*pow(timeStep,13)*pow(omega,13) + 1472*pow(timeStep,12)*pow(omega,12) + 8256*pow(timeStep,11)*pow(omega,11) + 11776*pow(timeStep,10)*pow(omega,10) + 39520*pow(timeStep,9)*pow(omega,9) + 46784*pow(timeStep,8)*pow(omega,8) + 95488*pow(timeStep,7)*pow(omega,7) + 87808*pow(timeStep,6)*pow(omega,6) + 105728*pow(timeStep,5)*pow(omega,5) + 68096*pow(timeStep,4)*pow(omega,4) + 44032*pow(timeStep,3)*pow(omega,3) + 16384*pow(timeStep,2)*pow(omega,2) + 4352*timeStep*omega + 512))/512,
        -(timeStep*omega*(pow(timeStep,15)*pow(omega,15) + 2*pow(timeStep,14)*pow(omega,14) + 42*pow(timeStep,13)*pow(omega,13) + 76*pow(timeStep,12)*pow(omega,12) + 660*pow(timeStep,11)*pow(omega,11) + 1048*pow(timeStep,10)*pow(omega,10) + 4840*pow(timeStep,9)*pow(omega,9) + 6448*pow(timeStep,8)*pow(omega,8) + 16944*pow(timeStep,7)*pow(omega,7) + 17696*pow(timeStep,6)*pow(omega,6) + 26208*pow(timeStep,5)*pow(omega,5) + 19264*pow(timeStep,4)*pow(omega,4) + 14784*pow(timeStep,3)*pow(omega,3) + 6272*pow(timeStep,2)*pow(omega,2) + 1920*timeStep*omega + 256))/256,
        -(timeStep*omega*(pow(timeStep,13)*pow(omega,13) + 2*pow(timeStep,12)*pow(omega,12) + 36*pow(timeStep,11)*pow(omega,11) + 64*pow(timeStep,10)*pow(omega,10) + 460*pow(timeStep,9)*pow(omega,9) + 696*pow(timeStep,8)*pow(omega,8) + 2528*pow(timeStep,7)*pow(omega,7) + 3040*pow(timeStep,6)*pow(omega,6) + 5808*pow(timeStep,5)*pow(omega,5) + 4960*pow(timeStep,4)*pow(omega,4) + 4672*pow(timeStep,3)*pow(omega,3) + 2304*pow(timeStep,2)*pow(omega,2) + 832*timeStep*omega + 128))/128,
        -(timeStep*omega*(pow(timeStep,11)*pow(omega,11) + 2*pow(timeStep,10)*pow(omega,10) + 30*pow(timeStep,9)*pow(omega,9) + 52*pow(timeStep,8)*pow(omega,8) + 296*pow(timeStep,7)*pow(omega,7) + 416*pow(timeStep,6)*pow(omega,6) + 1104*pow(timeStep,5)*pow(omega,5) + 1120*pow(timeStep,4)*pow(omega,4) + 1360*pow(timeStep,3)*pow(omega,3) + 800*pow(timeStep,2)*pow(omega,2) + 352*timeStep*omega + 64))/64,
        -(timeStep*omega*(pow(timeStep,9)*pow(omega,9) + 2*pow(timeStep,8)*pow(omega,8) + 24*pow(timeStep,7)*pow(omega,7) + 40*pow(timeStep,6)*pow(omega,6) + 168*pow(timeStep,5)*pow(omega,5) + 208*pow(timeStep,4)*pow(omega,4) + 352*pow(timeStep,3)*pow(omega,3) + 256*pow(timeStep,2)*pow(omega,2) + 144*timeStep*omega + 32))/32,
        -(timeStep*omega*(pow(timeStep,7)*pow(omega,7) + 2*pow(timeStep,6)*pow(omega,6) + 18*pow(timeStep,5)*pow(omega,5) + 28*pow(timeStep,4)*pow(omega,4) + 76*pow(timeStep,3)*pow(omega,3) + 72*pow(timeStep,2)*pow(omega,2) + 56*timeStep*omega + 16))/16,
        -(timeStep*omega*(pow(timeStep,5)*pow(omega,5) + 2*pow(timeStep,4)*pow(omega,4) + 12*pow(timeStep,3)*pow(omega,3) + 16*pow(timeStep,2)*pow(omega,2) + 20*timeStep*omega + 8))/8,
        -(timeStep*omega*(pow(timeStep,3)*pow(omega,3) + 2*pow(timeStep,2)*pow(omega,2) + 6*timeStep*omega + 4))/4,
        -(timeStep*omega*(timeStep*omega + 2))/2;

/*    Eigen::MatrixXd MAT_Pc_times_invM(2, 10);
    MAT_Pc_times_invM.row(0) = Pc_times_invM;
    MAT_Pc_times_invM.row(1) = Pc_times_invM;*/

    Eigen::MatrixXd bounds(10, 2);
    bounds.col(0) = 0.5*(Xz_max+Xz_min);
    bounds.col(1) = 0.5*(Yz_max+Yz_min);

    Eigen::Vector2d res = -Pc_times_invM.transpose()*bounds;

    return res;

  }

/*  void computeEndEffectorToWallConstraint(const Eigen::MatrixXd & M_ee,
                                          const Eigen::VectorXd & p_ee,
                                          const Eigen::VectorXd & n_obs,
                                          const Eigen::VectorXd & p_obs,
                                          Eigen::MatrixXd & A_constr,
                                          Eigen::VectorXd & b_constr);*/

  void computeInequalityConstraint(const Eigen::VectorXd & joint_pos, Eigen::MatrixXd & A_ineq, Eigen::VectorXd & b_ineq_L, Eigen::VectorXd & b_ineq_R);

  void computeMovingConstraint();

  // Takes a zmp_preview of dimension equal to prev
  Eigen::Vector2d computeXuStar(const Eigen::VectorXd & zmp_preview);

  // Takes a zmp_preview of dimension equal to (prev + N)
  Eigen::VectorXd computeXuStarPrediction(const Eigen::VectorXd & zmp_preview_full);

  void setParamsFromFile(const std::string & filepath);

  void setLogger() {

    logger_.setLoggingDirectory("./data");

    logger_.add(pred_zmp, "ZmpPredMPC", true, true);
    logger_.add(pred_left_pos, "LeftFootPredMPC", true, true);
    logger_.add(pred_right_pos, "RightFootPredMPC", true, true);
    logger_.add(com_pred_jacobians, "ComPosPredMPC", true, true);
    logger_.add(com_dot_pred_jacobians, "ComVelPredMPC", true, true);

//    logger_.add( new utils::Container<Eigen::VectorXd>(&pred_zmp,"zmp_pred"));
//    logger_.add( new utils::Container<Eigen::VectorXd>(&pred_left_pos,"left_foot_xy_pred"));
//    logger_.add( new utils::Container<Eigen::VectorXd>(&pred_right_pos,"right_foot_xy_pred"));
//
//    logger_.add( new utils::Container<Eigen::VectorXd>(&com_pred_jacobians,"com_pred_euler"));
//    logger_.add( new utils::Container<Eigen::VectorXd>(&xu_star,"tail"));

  //  logger_.add( new utils::Container<Eigen::VectorXd>(&com_dot_pred_jacobians,"com_vel_pred_euler"));

  }

};
