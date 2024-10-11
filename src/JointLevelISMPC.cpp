#include "JointLevelISMPC.hpp"

JointLevelISMPC::JointLevelISMPC(const int n_joints, const Eigen::MatrixXd & joint_lims, const Eigen::VectorXd & q_neutral_tot) {
  std::cout << "[JointLevelISMPC] initialization" << std::endl;

  q_neutral_ = q_neutral_tot.tail(n_joints);
  n_dof_ = n_joints + 6;
  nq_ = n_joints + 7;
  n_joints_ = n_joints;
  n_zmp_ = 2; // Dimension of ZMP position (x,y) and related CoM subspace

  this->setParamsFromFile("./mpc_params.toml");

  numVariables_ = n_dof_*N;
  numEqualityConstraints_ = n_zmp_ + 6*N;
  numInequalityConstraints_ = n_zmp_*N + n_joints_*n_joint_pos_constr_;// + N;

  A_stab_ = Eigen::MatrixXd::Zero(n_zmp_, numVariables_);
  b_stab_ = Eigen::VectorXd::Zero(n_zmp_);
  A_zmp_ = Eigen::MatrixXd::Zero(n_zmp_*(N), numVariables_);
  b_zmp_L_ = Eigen::VectorXd::Zero(A_zmp_.rows());
  b_zmp_R_ = Eigen::VectorXd::Zero(A_zmp_.rows());


  P_task_gen_ = Eigen::MatrixXd::Zero(N,N);
  p_task_gen_ = Eigen::MatrixXd::Zero(N,2);
  V_task_gen_ = Eigen::MatrixXd::Zero(N,N);
  v_task_gen_ = Eigen::MatrixXd::Zero(N,2);
  for (int row = 0; row<N; ++row) {
    for (int col = 0; col<=row; ++col) {
      P_task_gen_(row, col) = timeStep*timeStep/2.0 + (row-col)*timeStep*timeStep; //
      V_task_gen_(row, col) = timeStep;
    }
    p_task_gen_(row,0) = 1.0;
    p_task_gen_(row,1) = (1+row)*timeStep;
    v_task_gen_(row, 1) = 1.0;
  }


  kD_lim_L_ = -Eigen::VectorXd::Ones(b_zmp_L_.size())*footConstraintSquareWidth/2.0;
  kD_lim_R_ = Eigen::VectorXd::Ones(b_zmp_R_.size())*footConstraintSquareWidth/2.0;


  rk_ = Eigen::MatrixXd::Zero(n_zmp_*(N), n_zmp_);
  rk_prev_ = Eigen::MatrixXd::Zero(n_zmp_*(prev), n_zmp_);
  summation_ = Eigen::MatrixXd::Ones(n_zmp_, n_zmp_ *(N));
  summation_pred_ = Eigen::MatrixXd::Ones(n_zmp_, n_zmp_ *(prev));


  for(int i = 0; i<N; ++i)
    rk_.middleRows(i* n_zmp_, n_zmp_) = Eigen::MatrixXd::Identity(n_zmp_, n_zmp_);
  for(int i = 0; i<prev; ++i)
    rk_prev_.middleRows(i* n_zmp_, n_zmp_) = Eigen::MatrixXd::Identity(n_zmp_, n_zmp_);
  for(int i = 0; i<N; ++i)
    summation_.middleCols(i*n_zmp_,n_zmp_)= pow(exp(-omega*timeStep),i)*Eigen::MatrixXd::Identity(n_zmp_, n_zmp_);
  for(int i = 0; i<prev; ++i)
    summation_pred_.middleCols(i*n_zmp_,n_zmp_) = exp(-omega*timeStep*i)*Eigen::MatrixXd::Identity(n_zmp_, n_zmp_);

  selectX_ = Eigen::MatrixXd::Zero(N,2*N);
  selectY_ = Eigen::MatrixXd::Zero(N,2*N);

  for (int i = 0; i < N; ++i) {
    selectX_(i, i*2) = 1;
    selectY_(i, i*2+1) = 1;
  }

  // Joint position constraints
  A_joint_pos_ = Eigen::MatrixXd::Zero(n_joints*n_joint_pos_constr_, numVariables_);
/*  for(int r = 0; r<n_joint_pos_constr_; ++r)
    for(int c = 0; c<=r; ++c)
      A_joint_pos_.block(r*n_joints,c*n_dof_,n_joints,n_joints) = timeStep*Eigen::MatrixXd::Identity(n_joints,n_joints);*/
  b_joint_pos_L_ = utils::vrepeat(joint_lims.col(0), n_joint_pos_constr_);
  b_joint_pos_R_ = utils::vrepeat(joint_lims.col(1), n_joint_pos_constr_);

  current_vel_vector_ = Eigen::VectorXd::Zero(2*(N));

  // Compute moving constraint matrix that goes from 0 to 1 depending on the support foot. Zero for the right foot
  this->computeMovingConstraint();
  alpha_matrix_ = Eigen::MatrixXd::Zero(2*(N),2*(N));

  com_pred_jacobians = Eigen::VectorXd::Zero(3*N);

  COMStateToFeasibilityRegion_ = COMStateToFeasibilityRegion();

  this->setLogger();

  qp_solver_ptr_ = std::make_shared<labrob::qpsolvers::QPSolverEigenWrapper<double>>(
      std::make_shared<labrob::qpsolvers::HPIPMQPSolver>(
          numVariables_, numEqualityConstraints_, numInequalityConstraints_));

  // Initializations

  //------------------------------------------------------------------------------------------------------------------
  /* Compute matrix for regulation of joint positions
   * assuming [chest + minimal arms + legs] joint configuration */
  //------------------------------------------------------------------------------------------------------------------
/*  Eigen::MatrixXd joint_sel_matrix = Eigen::MatrixXd::Zero(n_joints, n_joints);
  joint_sel_matrix(0,0) = 20; // weight more the CHEST PITCH because moves too much otherwise
  joint_sel_matrix(1,1) = 1;
  joint_sel_matrix(2,2) = 1;
  joint_sel_matrix(3,3) = 10; // shoulder roll
  joint_sel_matrix(4,4) = 1;
  joint_sel_matrix(5,5) = 1;
  joint_sel_matrix(6,6) = 10; // shoulder roll
  joint_sel_matrix(7,7) = 1;
  joint_sel_matrix(8,8) = 1; // left hip yaw
  joint_sel_matrix(14,14) = 1; //right hip yaw

  Eigen::MatrixXd joint_sel_matrix_full = Eigen::MatrixXd::Zero(n_joints, n_dof_);
  joint_sel_matrix_full.rightCols(n_joints) = joint_sel_matrix;*/

  /* The whole matrix is the lower triangular repetition of joint_sel_matrix
   * times the JointLevelISMPC timeStep */
  joint_pos_matrix = Eigen::MatrixXd::Zero(N*n_joints, N*n_dof_);
  joint_pos_vector = Eigen::VectorXd(N*n_joints);
/*  for(int i = 0; i<N; ++i) {
    joint_pos_matrix.block(i*n_joints,i*n_dof_,(N-i)*n_joints,n_dof_) = timeStep*utils::vrepeat(joint_sel_matrix_full,N-i);
  }*/


  //------------------------------------------------------------------------------------------------------------------
  // Gains over prediction horizon into diagonal matrices
  //------------------------------------------------------------------------------------------------------------------
//  CoM_height_gains = reg_gains.CoM_height * Eigen::MatrixXd::Identity(N,N);
//  torso_orient_gains = utils::diagrepeat(reg_gains.torso_orient.asDiagonal().toDenseMatrix(),N);
//  swing_pose_gains = utils::diagrepeat(reg_gains.swing_pose.asDiagonal().toDenseMatrix(),N);

  // Kinematics
//  J_CoM = Eigen::MatrixXd::Zero(N,numVariables_);
//  z_CoM_pred_pos = Eigen::VectorXd::Zero(N);
  z_CoM_ref_pos = Eigen::VectorXd::Zero(N);
//  z_CoM_ref_vel = Eigen::VectorXd::Zero(N);

//  avrg_ang_velocity = Eigen::VectorXd::Zero(3*N);





}

void JointLevelISMPC::prepareMPC(const Data & data, const WalkState & walkstate) {

//  std::cout << "[JointLevelISMPC] Compute joint accelerations" << std::endl;
  costFunctionH = Eigen::MatrixXd::Zero(numVariables_, numVariables_);
  costFunctionF = Eigen::VectorXd::Zero(numVariables_);
  A_eq = Eigen::MatrixXd::Zero(numEqualityConstraints_, numVariables_);
  b_eq = Eigen::VectorXd::Zero(numEqualityConstraints_);
  A_ineq = Eigen::MatrixXd::Zero(numInequalityConstraints_,numVariables_);
  b_ineq_L = Eigen::VectorXd::Zero(numInequalityConstraints_);
  b_ineq_R = Eigen::VectorXd::Zero(numInequalityConstraints_);

  P_xyCoM = Eigen::MatrixXd::Zero(2*N,numVariables_);
  p_xyCoM = Eigen::VectorXd::Zero(2*N);

  V_xyCoM = Eigen::MatrixXd::Zero(2*N,numVariables_);
  v_xyCoM = Eigen::VectorXd::Zero(2*N);


  Eigen::Vector3d current_angmom = data.current_state.CoM.ang_vel;

//
//  angmom_in_zmp_(0) = current_angmom(1);
//  angmom_in_zmp_(1) = -current_angmom(0); // change sign because we need it for the ZMP
//  current_vel_vector_.head(2) = data.current_com_vel.head(2)/(omega*omega*timeStep);
//  current_vel_vector_.head(2) = data.current_com_vel.head(2)/(omega*omega*timeStep) + angmom_in_zmp_/(9.81*robot_mass*timeStep);


  //------------------------------------------------------------------------------------------------------------------
  // Build alpha matrix for the moving constraint, switching between left and right foot
  //------------------------------------------------------------------------------------------------------------------
  for(int i = 0; i<(N); ++i) {
    alpha_matrix_.block(i * 2, i * 2, 2, 2) = alpha_[walkstate.global_iter + (i) * timesRatio] * Eigen::Matrix2d::Identity();
  }

  //------------------------------------------------------------------------------------------------------------------
  // We have to extract the z part of the CoM kinematics and all the related quantities
  // Also set the average angular velocity and centroidal momentum matrix
  //------------------------------------------------------------------------------------------------------------------



  Eigen::MatrixXd P_zCoM = Eigen::MatrixXd::Zero(N,numVariables_);
  Eigen::VectorXd p_zCoM = Eigen::VectorXd::Zero(N);

  std::vector<Eigen::MatrixXd> zCoM_jacobians;
  std::vector<Eigen::MatrixXd> zCoM_jacobians_dot;
  for(int i = 0; i<N; ++i) {
    zCoM_jacobians.push_back(data.CoM_preds.jacobians[i].row(2));
    zCoM_jacobians_dot.push_back(data.CoM_preds.jacobians_dot[i].row(2));
  }
  Eigen::VectorXd zCoM_curr(1);
  zCoM_curr << data.current_state.CoM.pos(2);
  Eigen::VectorXd zCom_vel_curr(1);
  zCom_vel_curr << data.current_state.CoM.vel(2);
  buildTaskPosMatrix(zCoM_jacobians,zCoM_jacobians_dot, data.qdot, zCoM_curr , zCom_vel_curr, P_zCoM, p_zCoM);



  std::vector<Eigen::MatrixXd> xyCoM_jacobians;
  std::vector<Eigen::MatrixXd> xyCoM_jacobians_dot;

  for(int i = 0; i<N; ++i) {
    xyCoM_jacobians.push_back(data.CoM_preds.jacobians[i].topRows(2));
    xyCoM_jacobians_dot.push_back(data.CoM_preds.jacobians_dot[i].topRows(2));
  }
  buildTaskPosMatrix(xyCoM_jacobians, xyCoM_jacobians_dot, data.qdot, data.current_state.CoM.pos.head(2), data.current_state.CoM.vel.head(2), P_xyCoM, p_xyCoM);

//  V_xyCoM = Eigen::MatrixXd::Zero(2*N,numVariables_);
//  v_xyCoM = Eigen::VectorXd::Zero(2*N);
  buildTaskVelMatrix(xyCoM_jacobians, xyCoM_jacobians_dot, data.qdot, data.current_state.CoM.vel.head(2), V_xyCoM, v_xyCoM);

  Eigen::MatrixXd PN = (P_xyCoM.bottomRows(2) + 1/omega*V_xyCoM.bottomRows(2));
  Eigen::VectorXd pk = (p_xyCoM.tail(2) + 1/omega*v_xyCoM.tail(2));

  Eigen::MatrixXd P_feas_reg = COMStateToFeasibilityRegion_(0)*P_xyCoM + COMStateToFeasibilityRegion_(1)*V_xyCoM;
  Eigen::VectorXd p_feas_reg = COMStateToFeasibilityRegion_(0)*p_xyCoM + COMStateToFeasibilityRegion_(1)*v_xyCoM;


  Eigen::VectorXd feas_reg_ref = Eigen::VectorXd::Zero(2*N);
  Eigen::VectorXd zmp_ref_full = Eigen::VectorXd::Zero(n_zmp_*(data.zmp_refs.pos.size() + data.zmp_prev.size()));
  zmp_ref_full.head(n_zmp_*data.zmp_refs.pos.size()) = utils::vstack(data.zmp_refs.pos);
  zmp_ref_full.tail(n_zmp_*data.zmp_prev.size()) = utils::vstack(data.zmp_prev);
  for (int i = 0; i < N; ++i) {
    feas_reg_ref.segment(2*i,2) = computeXuStar(zmp_ref_full.segment(2*(N+i+1), 2*prev))
        + DesFeasibilityFromXzBounds(zmp_ref_full.segment(2*(i+1), 2*N));
  }

  for(int i = 0; i<N; ++i) {
//    J_CoM.block(i,i*n_dof_,1,n_dof_) = data.CoM_preds.jacobians[i].row(2);
//    z_CoM_pred_pos(i) = data.CoM_preds.pos[i](2);
    z_CoM_ref_pos(i) = data.CoM_refs.pos[i](2);
//    z_CoM_ref_vel(i) = data.CoM_refs.vel[i](2);
//    avrg_ang_velocity(3*i+2) = data.torso_orientation_refs.vel[i](2);
  }
  Eigen::MatrixXd Centroidal_Jacobian = utils::blkdiag(data.centroidal_momentum_preds.jacobians);

  Eigen::MatrixXd L_matrix = Eigen::MatrixXd::Zero(3*N, numVariables_);
  Eigen::VectorXd l_matrix = Eigen::VectorXd::Zero(3*N);
  buildTaskVelMatrix(data.centroidal_momentum_preds.jacobians,
                     data.centroidal_momentum_preds.jacobians_dot,
                     data.qdot,
                     current_angmom,
                     L_matrix,
                     l_matrix);


  //------------------------------------------------------------------------------------------------------------------
  // set the other task jacobians for torso orientation, swing foot and ZMP prediction matrix
  //------------------------------------------------------------------------------------------------------------------
  Eigen::MatrixXd P_torso_orient = Eigen::MatrixXd::Zero(3*N, numVariables_);
  Eigen::VectorXd p_torso_orient = Eigen::VectorXd::Zero(3*N);
  buildTaskPosMatrix(data.torso_orientation_preds.jacobians,
                     data.torso_orientation_preds.jacobians_dot,
                     data.qdot,
                     data.torso_orientation_preds.pos[0],
                     data.torso_orientation_preds.jacobians[0]*data.current_qdot,
                     P_torso_orient,
                     p_torso_orient);


  P_left_foot = Eigen::MatrixXd::Zero(6*N, numVariables_);
  p_left_foot = Eigen::VectorXd::Zero(6*N);
  buildTaskPosMatrix(data.left_foot_preds.jacobians,
                     data.left_foot_preds.jacobians_dot,
                     data.qdot,
                     data.left_foot_preds.pos[0],
                     data.left_foot_preds.jacobians[0]*data.current_qdot,
                     P_left_foot,
                     p_left_foot);
  P_right_foot = Eigen::MatrixXd::Zero(6*N, numVariables_);
  p_right_foot = Eigen::VectorXd::Zero(6*N);
  buildTaskPosMatrix(data.right_foot_preds.jacobians,
                     data.right_foot_preds.jacobians_dot,
                     data.qdot,
                     data.right_foot_preds.pos[0],
                     data.right_foot_preds.jacobians[0]*data.current_qdot,
                     P_right_foot,
                     p_right_foot);


  Eigen::MatrixXd V_left_foot = Eigen::MatrixXd::Zero(6*N, numVariables_);
  Eigen::VectorXd v_left_foot = Eigen::VectorXd::Zero(6*N);
  buildTaskVelMatrix(data.left_foot_preds.jacobians,
                     data.left_foot_preds.jacobians_dot,
                     data.qdot,
                     data.left_foot_preds.jacobians[0]*data.current_qdot,
                     V_left_foot,
                     v_left_foot);
  Eigen::MatrixXd V_right_foot = Eigen::MatrixXd::Zero(6*N, numVariables_);
  Eigen::VectorXd v_right_foot = Eigen::VectorXd::Zero(6*N);
  buildTaskVelMatrix(data.right_foot_preds.jacobians, data.right_foot_preds.jacobians_dot, data.qdot, data.right_foot_preds.jacobians[0]*data.current_qdot, V_right_foot,
                     v_right_foot);

  Eigen::MatrixXd P_right_hand = Eigen::MatrixXd::Zero(6*N, numVariables_);
  Eigen::VectorXd p_right_hand = Eigen::VectorXd::Zero(6*N);
  buildTaskPosMatrix(data.right_hand_preds.jacobians, data.right_hand_preds.jacobians_dot, data.qdot, data.right_hand_preds.pos[0], data.right_hand_preds.jacobians[0]*data.current_qdot, P_right_hand, p_right_hand);

  //Eigen::MatrixXd J_swing = utils::blkdiag(data.swing_foot_preds.jacobians);
  //Eigen::MatrixXd J_left_foot = utils::blkdiag(data.left_foot_preds.jacobians);
  //Eigen::MatrixXd J_right_foot = utils::blkdiag(data.right_foot_preds.jacobians);

  Eigen::MatrixXd L2ZMP = Eigen::MatrixXd::Zero(2,3);
  L2ZMP(0,1) = -1;
  L2ZMP(1,0) = 1;

  Eigen::MatrixXd Mk = this->joints2zmp(data.CoM_preds.jacobians, P_xyCoM);
//  Eigen::MatrixXd Mk_angmom = Mk;
  Eigen::MatrixXd Mk_angmom = this->joints2zmp_angular_momentum(data.CoM_preds.jacobians, data.centroidal_momentum_preds.jacobians, P_xyCoM);

  mk = Eigen::VectorXd::Zero(n_zmp_*N);
  mk << data.current_state.CoM.pos.head(2), p_xyCoM.head(n_zmp_*(N-1));
  mk += - utils::blkdiag(xyCoM_jacobians_dot)*utils::vstack(data.qdot)/(omega*omega);

  Eigen::VectorXd mk_angmom = mk + utils::diagrepeat(L2ZMP,N)*utils::blkdiag(data.centroidal_momentum_preds.jacobians_dot)*utils::vstack(data.qdot)/(9.81*robot_mass); // Add the Jdot*nu_k contribution

/*  Mk = Mk_angmom;
  mk = mk_angmom;*/

  //Eigen::MatrixXd Mk_angmom = this->joints2zmp_angular_momentum(data.CoM_preds.jacobians, data.CentroidalJacobian);

/*  std::cout << "Mk" << std::endl;
  std::cout << Mk << std::endl;
  std::cout << "mk" << std::endl;
  std::cout << mk.transpose() << std::endl;*/
  // ------------------------------------------------------------------------------------------------
  // Compute errors
  // ------------------------------------------------------------------------------------------------

  Eigen::VectorXd CoM_error = p_zCoM - z_CoM_ref_pos;
  Eigen::VectorXd zmp_error = mk_angmom - utils::vstack(data.zmp_refs.pos);
  Eigen::VectorXd torso_orient_error = Eigen::VectorXd::Zero(3*N);
  //Eigen::VectorXd swing_pose_error = utils::vstack(data.swing_foot_refs.pos) - utils::vstack(data.swing_foot_preds.pos);
  Eigen::VectorXd left_pose_error = p_left_foot - utils::vstack(data.left_foot_refs.pos);
  Eigen::VectorXd right_pose_error = p_right_foot - utils::vstack(data.right_foot_refs.pos);

  Eigen::VectorXd right_hand_pose_error = p_right_hand - utils::vstack(data.right_hand_refs.pos);

  // ensure proper orientation error
  for (int i = 0; i < N; ++i) {
    torso_orient_error.segment(i*3,3) = utils::angleSignedDistance(p_torso_orient.segment(i*3,3), data.torso_orientation_refs.pos[i]);
    //swing_pose_error.segment(i*6+3,3) = utils::angleSignedDistance(data.swing_foot_refs.pos[i].tail(3),data.swing_foot_preds.pos[i].tail(3));
    left_pose_error.segment(i*6+3,3) = utils::angleSignedDistance(p_left_foot.segment(i*6+3,3), data.left_foot_refs.pos[i].tail(3));
    right_pose_error.segment(i*6+3,3) = utils::angleSignedDistance(p_right_foot.segment(i*6+3,3), data.right_foot_refs.pos[i].tail(3));
    right_hand_pose_error.segment(i*6+3,3) = utils::angleSignedDistance(p_right_hand.segment(i*6+3,3), data.right_hand_refs.pos[i].tail(3));
  }


  // ------------------------------------------------------------------------------------------------
  // Assemble cost function terms and sum them with a regularization term
  // ------------------------------------------------------------------------------------------------

  // ------------------------------------------
  // JOINT ACCELERATIONS
  costFunctionH = qddot_weight * Eigen::MatrixXd::Identity(numVariables_, numVariables_);
  costFunctionF = Eigen::VectorXd::Zero(numVariables_);

  // ------------------------------------------
  // Feasibility region tracking
  costFunctionH += feas_weight * P_feas_reg.transpose() * P_feas_reg;
  costFunctionF += feas_weight * P_feas_reg.transpose() * (p_feas_reg - feas_reg_ref);

  // ------------------------------------------
  // COM HEIGHT
  costFunctionH += CoM_height_weight * P_zCoM.transpose() * P_zCoM;
  costFunctionF += CoM_height_weight * P_zCoM.transpose() * CoM_error;

  // ------------------------------------------
  // ZMP POSITION TRACKING


  // Fixed reference
  costFunctionH += zmp_weight * Mk_angmom.transpose() * Mk_angmom;
  costFunctionF += zmp_weight * Mk_angmom.transpose() * zmp_error;

  // ------------------------------------------
  // TORSO ORIENTATION TRACKING (FEEDFORWARD + POSITION ERROR)
  costFunctionH += torso_weight * P_torso_orient.transpose() * P_torso_orient;
  costFunctionF += torso_weight * P_torso_orient.transpose() * torso_orient_error;

  // ------------------------------------------
  // SWING FOOT TRACKING (FEEDFORWARD + POSITION ERROR)
  // Set to 0 xy tracking due to automatic footstep placement
  /*

  costFunctionH += swing_weight * J_swing.transpose()*utils::diagrepeat(auto_weights_swing,N)*J_swing;
  costFunctionF += - swing_weight * J_swing.transpose() *utils::diagrepeat(auto_weights_swing,N)* (utils::vstack(data.swing_foot_refs.vel) + swing_pose_gains*swing_pose_error);
  */

  std::vector<Eigen::Matrix6d> foot_weight_vec;

  for (int i = 0; i<N; ++i) {
    Eigen::Matrix6d foot_weight = Eigen::Matrix6d::Identity();
      /*if ((walkstate.mpc_iter + i + 1) % (S + D) >= S) { // I
        foot_weight(0, 0) = 0;
        foot_weight(1, 1) = 0;

      }*/
    foot_weight_vec.push_back(foot_weight);
  }

  costFunctionH += swing_weight * P_left_foot.transpose() * utils::blkdiag(foot_weight_vec) * P_left_foot;
  costFunctionF += swing_weight * P_left_foot.transpose() * utils::blkdiag(foot_weight_vec) * left_pose_error;
  
  costFunctionH += swing_weight * P_right_foot.transpose() * utils::blkdiag(foot_weight_vec)  * P_right_foot;
  costFunctionF += swing_weight * P_right_foot.transpose() * utils::blkdiag(foot_weight_vec)  * right_pose_error;

  Eigen::Matrix6d hand_weights_vec = Eigen::Vector6d(1,1,1,0,0,0).asDiagonal();

  // Cost function term for tracking of hand task
  // costFunctionH += 1 * P_right_hand.transpose() * utils::diagrepeat(hand_weights_vec,N) * P_right_hand;
  // costFunctionF += 1 * P_right_hand.transpose() * utils::diagrepeat(hand_weights_vec,N) * right_hand_pose_error;


  // ------------------------------------------
  // MINIMIZE CENTROIDAL ANGULAR MOMENTUM
  Eigen::Matrix3d Ldot_weight_matrix = Ldot_weight.asDiagonal();
  costFunctionH += Centroidal_Jacobian.transpose()*utils::diagrepeat(Ldot_weight_matrix,N)*Centroidal_Jacobian;
 // costFunctionF += Centroidal_Jacobian.transpose()*utils::diagrepeat(Ldot_weight_matrix,N)*utils::blkdiag(data.centroidal_momentum_preds.jacobians_dot) * utils::vrepeat(data.current_qdot, N);

  Eigen::Matrix3d avrg_ang_vel_weight_matrix = avrg_ang_vel_weight.asDiagonal();

  costFunctionH += L_matrix.transpose()*utils::diagrepeat(avrg_ang_vel_weight_matrix,N)*L_matrix;
  costFunctionF += L_matrix.transpose()*utils::diagrepeat(avrg_ang_vel_weight_matrix,N)*l_matrix;

  // ------------------------------------------
  // REGULATE SELECTED JOINT POSITION TO ZERO

  Eigen::MatrixXd joint_sel_matrix = Eigen::MatrixXd::Zero(n_joints_, n_joints_);
  joint_sel_matrix(0,0) = joint_regulation_weights_[0]; //20; // weight more the CHEST PITCH because moves too much otherwise
  joint_sel_matrix(1,1) = joint_regulation_weights_[1]; // CHEST YAW
  joint_sel_matrix(2,2) = joint_regulation_weights_[2]; // (L) SHOULDER PITCH
  joint_sel_matrix(3,3) = joint_regulation_weights_[3]; //10; // (L) SHOULDER ROLL
  joint_sel_matrix(4,4) = joint_regulation_weights_[4]; // (L) ELBOW PITCH
  joint_sel_matrix(5,5) = joint_regulation_weights_[5]; // (R) SHOULDER PITCH
  joint_sel_matrix(6,6) = joint_regulation_weights_[6]; //10; // (R) SHOULDER ROLL
  joint_sel_matrix(7,7) = joint_regulation_weights_[7]; // (R) ELBOW PITCH
  joint_sel_matrix(8,8) = joint_regulation_weights_[8]; // left hip yaw
  joint_sel_matrix(14,14) = joint_regulation_weights_[9]; //right hip yaw

  Eigen::MatrixXd joint_sel_matrix_full = Eigen::MatrixXd::Zero(n_joints_, n_dof_);
  joint_sel_matrix_full.rightCols(n_joints_) = Eigen::MatrixXd::Identity(n_joints_, n_joints_);

  std::vector<Eigen::MatrixXd> joint_jacobians;
  std::vector<Eigen::MatrixXd> joint_jacobians_dot;
  for(int i = 0; i<N; ++i) {
    joint_jacobians.push_back(joint_sel_matrix_full);
    joint_jacobians_dot.push_back(0*joint_sel_matrix_full);
  }
  joint_pos_vector.setZero();
  joint_pos_matrix.setZero();
  buildTaskPosMatrix(joint_jacobians,joint_jacobians_dot, data.qdot, data.current_q.tail(n_joints_), data.current_qdot.tail(n_joints_), joint_pos_matrix, joint_pos_vector);

  costFunctionH += joint_pos_weight * joint_pos_matrix.transpose() *utils::diagrepeat(joint_sel_matrix, N) * joint_pos_matrix;
  costFunctionF += joint_pos_weight * joint_pos_matrix.transpose() *utils::diagrepeat(joint_sel_matrix, N)* (joint_pos_vector - utils::vrepeat(q_neutral_,N));

  Eigen::MatrixXd joint_vel_matrix(joint_pos_matrix.rows(), joint_pos_matrix.cols());
  Eigen::VectorXd joint_vel_vector(joint_pos_vector.rows());

  buildTaskVelMatrix(joint_jacobians,joint_jacobians_dot, data.qdot, data.current_qdot.tail(n_joints_), joint_vel_matrix, joint_vel_vector);
  costFunctionH += joint_vel_weight * joint_vel_matrix.transpose() *utils::diagrepeat(joint_sel_matrix, N) * joint_vel_matrix;
  costFunctionF += joint_vel_weight * joint_vel_matrix.transpose() *utils::diagrepeat(joint_sel_matrix, N)* joint_vel_vector;

  A_joint_pos_ = joint_pos_matrix.topRows(n_joints_*n_joint_pos_constr_);

  // ------------------------------------------------------------------------------------------------
  // CONSTRAINTS
  // ------------------------------------------------------------------------------------------------

  // ------------------------------------------
  // input constraints (NOT USED)
  Eigen::MatrixXd jointLimConstrMatrix = Eigen::MatrixXd::Identity(numVariables_, numVariables_);
  Eigen::VectorXd jointLimUpperBound = 1 * Eigen::VectorXd::Ones(numVariables_);
  Eigen::VectorXd jointLimLowerBound = -1 * Eigen::VectorXd::Ones(numVariables_);
  // ------------------------------------------
  // Anticipative ZMP positions
  Eigen::VectorXd Xz_ant = utils::vstack(data.zmp_prev);
  // Compute stability and ZMP constraint

  this->computeStabilityConstraint(pk, PN, Xz_ant.head(n_zmp_*prev));
  this->computeZMPConstraint(mk_angmom, Mk_angmom, data.zmp_refs.pos, data);
  // this->computeZMPConstraint(mk, Mk, data.zmp_refs.pos, data);
//  this->computeAutoZMPConstraint(data.current_com_pos.head(2),right_foot_pos_k,left_foot_pos_k, Mk, M_right_foot,M_left_foot, data);


  // ------------------------------------------
  // Inequality constraints
  // Allocate the final data structures
  
  // fill the data structures with joint limits and ZMP constraint
  this->computeInequalityConstraint(data.current_q.tail(n_joints_), A_ineq, b_ineq_L, b_ineq_R);
  // ------------------------------------------

/*  A_ineq*=0;
  b_ineq_L*=0;
  b_ineq_R*=0;*/


  Eigen::MatrixXd A_support_foot = Eigen::MatrixXd::Zero(6*N, numVariables_);
  Eigen::VectorXd b_support_foot = Eigen::VectorXd::Zero(6*N);

  Eigen::MatrixXd A_double_support = Eigen::MatrixXd::Zero(2*N, numVariables_);
  Eigen::VectorXd b_double_support = Eigen::VectorXd::Zero(2*N);


  for (int i = 0; i < N; ++i) {
    if (walkstate.predicted_support_foot.at(i+1) == Foot::RIGHT) {
      A_support_foot.middleRows(6*i, 6) = V_right_foot.middleRows(6*i, 6);
      b_support_foot.segment(i*6,6) = - v_right_foot.middleRows(6*i, 6);

      if((walkstate.mpc_iter+i+1)%(S+D)>=S) {
        A_double_support.middleRows(i*2, 2) = V_left_foot.middleRows(6*i, 2);
        b_double_support.segment(i*2,2) =  - v_left_foot.middleRows(6*i, 2);
      }

    } else {
      A_support_foot.middleRows(6*i, 6) = V_left_foot.middleRows(6*i, 6);
      b_support_foot.segment(i*6,6) = - v_left_foot.middleRows(6*i, 6);

      if((walkstate.mpc_iter+i+1)%(S+D)>=S) {
        A_double_support.middleRows(i*2, 2) = V_right_foot.middleRows(6*i, 2);
        b_double_support.segment(i*2,2) =  - v_right_foot.middleRows(6*i, 2);
      }
    }
  }

  costFunctionH += double_support_weight*A_double_support.transpose()*A_double_support;
  costFunctionF += double_support_weight*A_double_support.transpose()*(-b_double_support);

  A_eq << A_stab_, A_support_foot;
  b_eq << b_stab_, b_support_foot;


  // For logging
  Mk_log_ = Mk_angmom;
  mk_log_ = mk_angmom;
  P_xyCoM_log_ = P_xyCoM;
  p_xyCoM_log_ = p_xyCoM;


}  

Eigen::MatrixXd JointLevelISMPC::computeQAcc(const Data & data, const WalkState & walkstate) {

  // ------------------------------------------------------------------------------------------------
  // SOLVE QP, GET SOLUTION, LOGGING
  // ------------------------------------------------------------------------------------------------
  qp_solver_ptr_->solve(
      costFunctionH,
      costFunctionF,
      A_eq,
      b_eq,
      A_ineq,
      b_ineq_L,
      b_ineq_R);

  Eigen::VectorXd qAcc_all = qp_solver_ptr_->get_solution();

  double tolerance = 1e-3;
  bool eq_constr_satisfaction = (A_eq * qAcc_all - b_eq).lpNorm<Eigen::Infinity>() < tolerance;
  Eigen::VectorXd leftRes = b_ineq_L - A_ineq * qAcc_all;
  Eigen::VectorXd rightRes = A_ineq * qAcc_all - b_ineq_R;
  bool ineq_constr_satisfaction = true;
  for (int i = 0; i < leftRes.rows(); ++i) {
    if (leftRes(i) > tolerance || rightRes(i) > tolerance) ineq_constr_satisfaction = false;
  }
  if (!eq_constr_satisfaction || !ineq_constr_satisfaction) {
    simulation_outcome = 1;
    std::cout << "Infeasibility!!" << std::endl;
  }

  Eigen::MatrixXd QAcc(n_dof_, N);
  for (int i = 0; i < N; i++) QAcc.col(i) = qAcc_all.segment(i*n_dof_, n_dof_);

  // LOG
  //current_vel_vector_.head(2) +=  angmom_in_zmp_/(9.81*robot_mass*timeStep);

  pred_zmp = Mk_log_*qAcc_all + mk_log_;
  pred_left_pos = P_left_foot*qAcc_all + p_left_foot;
  pred_right_pos = P_right_foot*qAcc_all + p_right_foot;

  com_pred_jacobians = P_xyCoM_log_*qAcc_all + p_xyCoM_log_;

  com_dot_pred_jacobians = V_xyCoM*qAcc_all + v_xyCoM;


  logger_.logList();

  return QAcc;
}

void JointLevelISMPC::buildTaskPosMatrix(const std::vector<Eigen::MatrixXd> & jacobians,
                                         const std::vector<Eigen::MatrixXd> & jacobians_dot,
                                         const std::vector<Eigen::VectorXd> & qdot,
                                         const Eigen::VectorXd & task_pos,
                                         const Eigen::VectorXd & task_vel,
                                         Eigen::MatrixXd & P_task,
                                         Eigen::VectorXd & p_task) {
  // Find out task dimension
  int tdim = task_vel.size();

  for (int row = 0; row<N; ++row) {
    for (int col = 0; col<=row; ++col) {
      P_task.block(row*tdim, col*n_dof_, tdim, n_dof_) = P_task_gen_(row, col) * jacobians[col];
      p_task.segment(row*tdim, tdim) += P_task_gen_(row, col)*jacobians_dot[col]*qdot[col];
    }

    p_task.segment(row*tdim, tdim) += p_task_gen_(row,0)*task_pos + p_task_gen_(row,1)*task_vel;
  }

}

void JointLevelISMPC::buildTaskVelMatrix(const std::vector<Eigen::MatrixXd> & jacobians,
                                         const std::vector<Eigen::MatrixXd> & jacobians_dot,
                                         const std::vector<Eigen::VectorXd> & qdot,
                                         const Eigen::VectorXd & task_vel,
                                         Eigen::MatrixXd & V_task,
                                         Eigen::VectorXd & v_task) {
  // Find out task dimension
  int tdim = task_vel.size();

  for (int row = 0; row<N; ++row) {
    for (int col = 0; col<=row; ++col) {
      V_task.block(row*tdim, col*n_dof_, tdim, n_dof_) = V_task_gen_(row, col) * jacobians[col];
      v_task.segment(row*tdim, tdim) += V_task_gen_(row, col)*jacobians_dot[col]*qdot[col];
    }
    v_task.segment(row*tdim, tdim) += v_task_gen_(row,1)*task_vel;
  }

}


void JointLevelISMPC::computeZMPConstraint(const Eigen::VectorXd & mk,
                                           const Eigen::MatrixXd & Mk,
                                           const std::vector<Eigen::Vector2d> & zmp_pos,
                                           const Data & data) {

  Eigen::VectorXd F_kp1 = utils::vstack(zmp_pos);

  A_zmp_ = Mk;
  b_zmp_L_ = kD_lim_L_ + F_kp1 - mk;
  b_zmp_R_ = kD_lim_R_ + F_kp1 - mk;
}

void JointLevelISMPC::computeAutoZMPConstraint(const Eigen::VectorXd & xck,
                                               const Eigen::VectorXd & xrk,
                                               const Eigen::VectorXd & xlk,
                                               const Eigen::MatrixXd & Mk,
                                               const Eigen::MatrixXd & Mr,
                                               const Eigen::MatrixXd & Ml,
                                               const Data & data) {

  Eigen::VectorXd current_vel_vector = Eigen::VectorXd::Zero(2*(N));
 // current_vel_vector.head(2) = data.current_com_vel.head(2)/(omega*omega*timeStep);

  A_zmp_ = Mk - alpha_matrix_*Mr - (Eigen::MatrixXd::Identity(2*(N), 2*(N))-alpha_matrix_)*Ml;
  b_zmp_L_ = kD_lim_L_ + alpha_matrix_*xrk + (Eigen::MatrixXd::Identity(2*(N), 2*(N))-alpha_matrix_)*xlk - rk_*xck - current_vel_vector;// - disturbance_over_zmp;
  b_zmp_R_ = kD_lim_R_ + alpha_matrix_*xrk + (Eigen::MatrixXd::Identity(2*(N), 2*(N))-alpha_matrix_)*xlk - rk_*xck - current_vel_vector;// - disturbance_over_zmp;

}

void JointLevelISMPC::computeStabilityConstraint(const Eigen::VectorXd & pk,
                                                 const Eigen::MatrixXd & PN,
                                                 const Eigen::VectorXd & Xz_ant) {
  xu_star = computeXuStar(Xz_ant);
  A_stab_ = PN;
  b_stab_ = xu_star - pk;
}

Eigen::Vector2d JointLevelISMPC::computeXuStar(const Eigen::VectorXd & zmp_preview) {
  return (1 - exp(-omega*timeStep))*summation_pred_*zmp_preview + exp(-omega*timeStep*prev)*zmp_preview.tail(n_zmp_);
}

Eigen::VectorXd JointLevelISMPC::computeXuStarPrediction(const Eigen::VectorXd & zmp_preview_full) {
  Eigen::VectorXd xu_star_preview = Eigen::VectorXd::Zero(2*N);
  for(int i=0; i<N; ++i){
    xu_star_preview.segment(2*i, 2) = computeXuStar(zmp_preview_full.segment(n_zmp_*i, n_zmp_*prev));
  }
  return xu_star_preview;
}


Eigen::MatrixXd JointLevelISMPC::joints2zmp(const std::vector<Eigen::MatrixXd> & CoM_jacobians, const Eigen::MatrixXd & P_xyCoM) const {


  Eigen::MatrixXd Mk = Eigen::MatrixXd::Zero(n_zmp_*N, numVariables_);
  for(int r = 0; r<N; ++r) {
    Mk.block(r * n_zmp_, r * n_dof_, n_zmp_, n_dof_) = - (CoM_jacobians[r].topRows(n_zmp_))/(omega*omega);
  }
  // Insert the CoM position component
  Mk.bottomLeftCorner((N-1)*n_zmp_, (N-1)*n_dof_) += P_xyCoM.topLeftCorner((N-1)*n_zmp_, (N-1)*n_dof_);

  return Mk;
}

Eigen::MatrixXd JointLevelISMPC::joints2zmp_angular_momentum(const std::vector<Eigen::MatrixXd> & CoM_jacobians,
                                                             const std::vector<Eigen::MatrixXd> & CMM,
                                                             const Eigen::MatrixXd & P_xyCoM) const  {
  Eigen::MatrixXd Mk = Eigen::MatrixXd::Zero(n_zmp_*N, numVariables_);
  for(int r = 0; r<N; ++r) {
    Eigen::MatrixXd Ak = Eigen::MatrixXd::Zero(2,n_dof_);
    Ak.row(0) = -CMM[r].row(1);
    Ak.row(1) = CMM[r].row(0);
    Mk.block(r * n_zmp_, r * n_dof_, n_zmp_, n_dof_) = - (CoM_jacobians[r].topRows(n_zmp_))/(omega*omega) + Ak/(robot_mass*9.81) ;
  }
  // Insert the CoM position component
  Mk.bottomLeftCorner((N-1)*n_zmp_, (N-1)*n_dof_) += P_xyCoM.topLeftCorner((N-1)*n_zmp_, (N-1)*n_dof_);

  return Mk;

}

void JointLevelISMPC::computeInequalityConstraint(const Eigen::VectorXd & joint_pos, Eigen::MatrixXd & A_ineq, Eigen::VectorXd & b_ineq_L, Eigen::VectorXd & b_ineq_R) {

  if (n_joint_pos_constr_ == 0) {
    A_ineq.middleRows(0, A_zmp_.rows()) = A_zmp_;
    b_ineq_L.segment(0, b_zmp_L_.size()) = b_zmp_L_;;
    b_ineq_R.segment(0, b_zmp_R_.size()) = b_zmp_R_;
  }
  else {
    A_ineq.middleRows(0, A_zmp_.rows()) = A_zmp_;
    A_ineq.middleRows(A_zmp_.rows(), A_joint_pos_.rows()) = A_joint_pos_;

    Eigen::VectorXd qk_stack = joint_pos_vector.head(n_joints_*n_joint_pos_constr_);

    b_ineq_L.segment(0, b_zmp_L_.size()) = b_zmp_L_;
    b_ineq_L.segment(b_zmp_L_.size(), b_joint_pos_L_.size()) = (b_joint_pos_L_ - qk_stack);
    b_ineq_R.segment(0, b_zmp_R_.size()) = b_zmp_R_;
    b_ineq_R.segment(b_zmp_R_.size(), b_joint_pos_R_.size()) = (b_joint_pos_R_ - qk_stack);
  }
}

void JointLevelISMPC::setParamsFromFile(const std::string & filepath) {

  toml::table config = toml::parse_file(filepath);

  n_joint_pos_constr_ = config["num_pos_constr"].as_integer()->get();

  // Cost function weights
  avrg_ang_vel_weight(0) = config.at_path("weights.avrg_ang_vel[0]").as_floating_point()->get();
  avrg_ang_vel_weight(1) = config.at_path("weights.avrg_ang_vel[1]").as_floating_point()->get();
  avrg_ang_vel_weight(2) = config.at_path("weights.avrg_ang_vel[2]").as_floating_point()->get();

  Ldot_weight(0) = config.at_path("weights.L_dot[0]").as_floating_point()->get();
  Ldot_weight(1) = config.at_path("weights.L_dot[1]").as_floating_point()->get();
  Ldot_weight(2) = config.at_path("weights.L_dot[2]").as_floating_point()->get();

  joint_pos_weight = config.at_path("weights.joint_pos").as_floating_point()->get();
  torso_weight = config.at_path("weights.torso_orient").as_floating_point()->get();
  swing_weight = config.at_path("weights.swing_pose").as_floating_point()->get();
  zmp_weight = config.at_path("weights.zmp_pos").as_floating_point()->get();
  feas_weight = config.at_path("weights.feas_pos").as_floating_point()->get();
  CoM_height_weight = config.at_path("weights.CoM_height").as_floating_point()->get();
  footstep_weight = config.at_path("weights.footstep").as_floating_point()->get();
  double_support_weight = config.at_path("weights.double_support").as_floating_point()->get();
  joint_vel_weight = config.at_path("weights.qdot").as_floating_point()->get();
  qddot_weight = config.at_path("weights.qddot").as_floating_point()->get();
  

  // Regulation gains for the feedback term
/*  reg_gains.CoM_height = config.at_path("gains.CoM_height").as_floating_point()->get();
  reg_gains.torso_orient(0) = config.at_path("gains.torso_orient[0]").as_floating_point()->get();
  reg_gains.torso_orient(1) = config.at_path("gains.torso_orient[1]").as_floating_point()->get();
  reg_gains.torso_orient(2) = config.at_path("gains.torso_orient[2]").as_floating_point()->get();

  reg_gains.swing_pose(0) = config.at_path("gains.swing_position[0]").as_floating_point()->get();
  reg_gains.swing_pose(1) = config.at_path("gains.swing_position[1]").as_floating_point()->get();
  reg_gains.swing_pose(2) = config.at_path("gains.swing_position[2]").as_floating_point()->get();

  reg_gains.swing_pose(3) = config.at_path("gains.swing_orient[0]").as_floating_point()->get();
  reg_gains.swing_pose(4) = config.at_path("gains.swing_orient[1]").as_floating_point()->get();
  reg_gains.swing_pose(5) = config.at_path("gains.swing_orient[2]").as_floating_point()->get();*/

  for(int i = 0; i<10; ++i){
    joint_regulation_weights_.push_back(config.at_path("weights.joint_regulation["+std::to_string(i)+"]").as_integer()->get());
  }

}

void JointLevelISMPC::computeMovingConstraint() {
  double alpha_start = 0.5;
  for(int step = 0; step<50; ++step) {
    double alpha_end = step % 2;

    for (int k = 0; k < singleSupportSamples; ++k) {
      alpha_.emplace_back(alpha_start);
    }
    for (int k = 0; k < doubleSupportSamples; ++k) {
      double tau = (k/(double)doubleSupportSamples);

      double alpha_k = alpha_start + tau*(alpha_end-alpha_start);
      alpha_.emplace_back(alpha_k);
    }
    alpha_start = alpha_end;
  }
}


  std::tuple<double, double, double, double> JointLevelISMPC::getFeasibilityRegionBounds(const Data & data) {

    Eigen::VectorXd Xc_ddot_min, Xc_ddot_max, Yc_ddot_min, Yc_ddot_max;

    Eigen::MatrixXd select_x = Eigen::MatrixXd::Zero(N, N*n_zmp_);
    Eigen::MatrixXd select_y = Eigen::MatrixXd::Zero(N, N*n_zmp_);
    for(int i = 0; i<N; ++i) {
      select_x(i, 2*i) = 1;
      select_y(i, 2*i+1) = 1;
    }

    Eigen::MatrixXd Mk_ct = Eigen::MatrixXd::Zero(N, N);

    for(int r = 0; r<N; ++r) {
      Mk_ct(r, r) = - 1/(omega*omega);
    }
    // Insert the CoM position component
    Mk_ct.bottomLeftCorner(N-1, N-1) += P_task_gen_.topLeftCorner((N-1), (N-1));

    Eigen::MatrixXd PN_ct = (P_task_gen_.bottomRows(1) + 1/omega*V_task_gen_.bottomRows(1));

    Xc_ddot_min = Mk_ct.inverse() * select_x * (kD_lim_L_ + utils::vstack(data.zmp_refs.pos) - mk);
    Xc_ddot_max = Mk_ct.inverse() * select_x * (kD_lim_R_ + utils::vstack(data.zmp_refs.pos) - mk);
    Yc_ddot_min = Mk_ct.inverse() * select_y * (kD_lim_L_ + utils::vstack(data.zmp_refs.pos) - mk);
    Yc_ddot_max = Mk_ct.inverse() * select_y * (kD_lim_R_ + utils::vstack(data.zmp_refs.pos) - mk);

    return std::make_tuple(
      - (PN_ct * Xc_ddot_min)(0) + xu_star(0),
      - (PN_ct * Xc_ddot_max)(0) + xu_star(0),
      - (PN_ct * Yc_ddot_min)(0) + xu_star(1),
      - (PN_ct * Yc_ddot_max)(0) + xu_star(1));
  }

