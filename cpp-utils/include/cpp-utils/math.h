#pragma once

#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <vector>


namespace utils {

    inline bool isEven(int n) { return not (n % 2); }

    inline bool isOdd(int n) { return (n % 2); }

    inline double wrapToPi(double angle){
      double ret=angle;
      while(ret>M_PI)
        ret-=2*M_PI;

      while(ret<=-M_PI)
        ret+=2*M_PI;

      return ret;
    }

    inline Eigen::MatrixXd matrixPower(const Eigen::MatrixXd& A, int exp){
      Eigen::MatrixXd result = Eigen::MatrixXd::Identity(A.rows(),A.cols());

      for (int i=0; i<exp;++i)
        result *= A;
      return result;
    }

    inline double sign(double x){
      if(x>0) return +1;
      if(x<0) return -1;
      return -1;
    }

// Elementary rotation around xyz axes
    inline Eigen::Matrix3d rotx(double ax){
      Eigen::Matrix3d rx = Eigen::Matrix3d::Zero();

      rx(0,0) = 1;

      rx(1,1) = cos(ax);
      rx(1,2) = -sin(ax);

      rx(2,1) = sin(ax);
      rx(2,2) = cos(ax);

      return rx;
    }

    inline Eigen::Matrix3d roty(double ay){
      Eigen::Matrix3d ry = Eigen::Matrix3d::Zero();

      ry(0,0) = cos(ay);
      ry(0,2) = sin(ay);

      ry(1,1) = 1;

      ry(2,0) = -sin(ay);
      ry(2,2) = cos(ay);

      return ry;
    }

    inline Eigen::Matrix3d rotz(double az){
      Eigen::Matrix3d rz = Eigen::Matrix3d::Zero();

      rz(0,0) = cos(az);
      rz(0,1) = -sin(az);

      rz(1,0) = sin(az);
      rz(1,1) = cos(az);

      rz(2,2) = 1;
      return rz;
    }

    inline double angleSignedDistance(double a, double b){
      double d = a - b;
      while(d >  M_PI) d = d - 2.0*M_PI;
      while(d < -M_PI) d = d + 2.0*M_PI;

      return d;
    }

    inline Eigen::Matrix3d rot(const Eigen::Vector3d & eul){

      Eigen::Matrix3d rx = rotx(eul(0));
      Eigen::Matrix3d ry = roty(eul(1));
      Eigen::Matrix3d rz = rotz(eul(2));

      return rz*ry*rx;
    }

    inline Eigen::Matrix2d rot(const double angle){

      Eigen::Matrix2d R;
      const double c = cos(angle);
      const double s = sin(angle);
      R(0,0) = c;
      R(0,1) = -s;
      R(1,0) = s;
      R(1,1) = c;

      return R;
    }

    inline Eigen::Vector3d angleSignedDistance(const Eigen::Vector3d & a, const Eigen::Vector3d & b) {
      Eigen::Matrix3d Ra = rot(a);
      Eigen::Matrix3d Rb = rot(b);

      Eigen::Matrix3d Rdiff = Rb.transpose() * Ra;
      auto aa = Eigen::AngleAxisd(Rdiff);
      return aa.angle() * Ra * aa.axis();
    }

    inline Eigen::Vector3d getRPY(const Eigen::Matrix3d & rotMatrix) {
      Eigen::Vector3d RPY;
      RPY << atan2(rotMatrix(2, 1), rotMatrix(2, 2)),
          atan2(-rotMatrix(2, 0), sqrt(rotMatrix(2, 1) * rotMatrix(2, 1) + rotMatrix(2, 2) * rotMatrix(2, 2))),
          atan2(rotMatrix(1, 0), rotMatrix(0, 0));
      return RPY;
    }



    inline Eigen::Matrix4d v2t(const Eigen::VectorXd & v){

      Eigen::Matrix4d m = Eigen::Matrix4d::Identity();

      Eigen::Vector3d eul = v.head(3);
      Eigen::Matrix3d r = utils::rot(eul);

      m.block<3,3>(0,0) = r;
      m.block<3,1>(0,3) = v.tail(3);

      return m;
    }

    inline Eigen::VectorXd t2v(const Eigen::Matrix4d & m) {
      Eigen::VectorXd v(6);

      double beta = atan2( m(0,2), sqrt(pow(m(0,0), 2) + pow(m(0,1), 2)) );
      double alpha = atan2( -m(1,2)/cos(beta), m(2,2)/cos(beta) );
      double gamma = atan2( -m(0,1)/cos(beta), m(0,0)/cos(beta) );

      v(0) = alpha;
      v(1) = beta;
      v(2) = gamma;
      v(3) = m(0,3);
      v(4) = m(1,3);
      v(5) = m(2,3);

      return v;
    }

// Express v2 in the frame of v1
    inline Eigen::VectorXd vvRel(const Eigen::VectorXd & v2, const Eigen::VectorXd & v1) {
      return utils::t2v(utils::v2t(v1).inverse() * utils::v2t(v2));
    }

// Check if an element is in a vector
    template<typename T>
    bool is_in_vector(const std::vector<T> & vector, const T & elt) {
      return vector.end() != std::find(vector.begin(),vector.end(),elt);
    }

// arrange elements in a block diagonal matrix
    template<typename MatrixEigen>
    Eigen::MatrixXd blkdiag(const std::vector<MatrixEigen> & a) {
      int count = a.size();
      int block_rows = a[0].rows();
      int block_cols = a[0].cols();
      Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(block_rows * count, block_cols * count);
      for (int i = 0; i < count; ++i)
        bdm.block(i * block_rows, i * block_cols, block_rows, block_cols) = a[i];

      return bdm;
    }

// Repeat in a vertical stack for a certain amount of times
    template<typename Derived>
    Eigen::MatrixXd vrepeat(const Eigen::MatrixBase<Derived> & a, int num) {
      int block_rows = a.rows();
      Eigen::MatrixXd vstack = Eigen::MatrixXd::Zero(block_rows * num, a.cols());
      for (int i = 0; i < num; ++i)
        vstack.middleRows(i * block_rows, block_rows) = a;
      return vstack;
    }

    inline Eigen::MatrixXd vrepeat(const Eigen::MatrixXd & a, int num) {
      int block_rows = a.rows();
      Eigen::MatrixXd vstack = Eigen::MatrixXd::Zero(block_rows * num, a.cols());
      for (int i = 0; i < num; ++i)
        vstack.middleRows(i * block_rows, block_rows) = a;
      return vstack;
    }


// Repeat elements diagonally for a certain amount of times
    template<typename Derived>
    Eigen::MatrixXd diagrepeat(const Eigen::DenseBase<Derived> & a, int num) {
      int block_rows = a.rows();
      int block_cols = a.cols();
      Eigen::MatrixXd diag = Eigen::MatrixXd::Zero(num*block_rows, num*block_cols);
      for (int i = 0; i < num; ++i)
        diag.block(i * block_rows,i * block_cols, block_rows, block_cols) = a;
      return diag;
    }

// Stack elements vertically from std::vector
    template<typename VectorEigen>
    Eigen::VectorXd vstack(const std::vector<VectorEigen> & a) {
      int num = a.size();
      int block_rows = a[0].rows();
      Eigen::VectorXd vstack = Eigen::VectorXd::Zero(block_rows * num);
      for (int i = 0; i < num; ++i)
        vstack.segment(i * block_rows, block_rows) = a[i];
      return vstack;
    }

    inline Eigen::MatrixXd vstack(const std::vector<Eigen::MatrixXd> & a) {
      int num = a.size();
      int block_rows = a[0].rows();
      Eigen::MatrixXd vstack = Eigen::MatrixXd::Zero(block_rows * num, a[0].cols());
      for (int i = 0; i < num; ++i)
        vstack.middleRows(i * block_rows, block_rows) = a[i];
      return vstack;
    }

    static Eigen::MatrixXd quinticPolynomialEE(int N, // number of samples
                                               double tf, // final time
                                               const Eigen::VectorXd& x0, // initial state
                                               const Eigen::VectorXd& xf, // final state
                                               const Eigen::VectorXd& dx0, // initial velocity
                                               const Eigen::VectorXd& dxf, // final velocity
                                               const Eigen::VectorXd& ddx0, // initial acceleration
                                               const Eigen::VectorXd& ddxf) { // final acceleration


      int n_dim = x0.size();

      Eigen::MatrixXd poly_interp = Eigen::MatrixXd::Zero(3*n_dim, N);

      Eigen::VectorXd boundary_conditions_stack(6*n_dim);
      boundary_conditions_stack << x0, dx0, ddx0, xf, dxf, ddxf;

      Eigen::MatrixXd constr_matrix = Eigen::MatrixXd::Zero(6*n_dim, 6*n_dim);

      Eigen::MatrixXd I_n = Eigen::MatrixXd::Identity(n_dim, n_dim);

      constr_matrix.block(0,0,3*n_dim,3*n_dim) = Eigen::MatrixXd::Identity(3*n_dim, 3*n_dim);

      constr_matrix.block(3*n_dim, 0, n_dim, n_dim) = I_n;
      constr_matrix.block(3*n_dim, n_dim, n_dim, n_dim) = tf * I_n;
      constr_matrix.block(3*n_dim, 2*n_dim, n_dim, n_dim) = pow(tf, 2) * I_n;
      constr_matrix.block(3*n_dim, 3*n_dim, n_dim, n_dim) = pow(tf, 3) * I_n;
      constr_matrix.block(3*n_dim, 4*n_dim, n_dim, n_dim) = pow(tf, 4) * I_n;
      constr_matrix.block(3*n_dim, 5*n_dim, n_dim, n_dim) = pow(tf, 5) * I_n;

      constr_matrix.block(4*n_dim, n_dim, n_dim, n_dim) = I_n;
      constr_matrix.block(4*n_dim, 2*n_dim, n_dim, n_dim) = 2 * tf * I_n;
      constr_matrix.block(4*n_dim, 3*n_dim, n_dim, n_dim) = 3 * pow(tf, 2) * I_n;
      constr_matrix.block(4*n_dim, 4*n_dim, n_dim, n_dim) = 4 * pow(tf, 3) * I_n;
      constr_matrix.block(4*n_dim, 5*n_dim, n_dim, n_dim) = 5 * pow(tf, 4) * I_n;

      constr_matrix.block(5*n_dim, 2*n_dim, n_dim, n_dim) = 2 * I_n;
      constr_matrix.block(5*n_dim, 3*n_dim, n_dim, n_dim) = 6 * tf * I_n;
      constr_matrix.block(5*n_dim, 4*n_dim, n_dim, n_dim) = 12 * pow(tf, 2) * I_n;
      constr_matrix.block(5*n_dim, 5*n_dim, n_dim, n_dim) = 20 * pow(tf, 3) * I_n;

      Eigen::VectorXd coeffs = constr_matrix.inverse() * boundary_conditions_stack;

      for(int i=0; i<N; ++i) {
        double t = tf/(N-1) * i;
        poly_interp.col(i).head(n_dim) = coeffs.segment(0, n_dim) +
                                         coeffs.segment(n_dim, n_dim) * t +
                                         coeffs.segment(2*n_dim, n_dim) * pow(t, 2) +
                                         coeffs.segment(3*n_dim, n_dim) * pow(t, 3) +
                                         coeffs.segment(4*n_dim, n_dim) * pow(t, 4) +
                                         coeffs.segment(5*n_dim, n_dim) * pow(t, 5);

        poly_interp.col(i).segment(n_dim,n_dim) =
            coeffs.segment(n_dim, n_dim) +
            coeffs.segment(2*n_dim, n_dim) * t * 2 +
            coeffs.segment(3*n_dim, n_dim) * pow(t, 2) * 3 +
            coeffs.segment(4*n_dim, n_dim) * pow(t, 3) * 4 +
            coeffs.segment(5*n_dim, n_dim) * pow(t, 4) * 5;

        poly_interp.col(i).tail(n_dim) =
            coeffs.segment(n_dim, n_dim) +
            coeffs.segment(2*n_dim, n_dim) * 2 +
            coeffs.segment(3*n_dim, n_dim) * t * 6 +
            coeffs.segment(4*n_dim, n_dim) * pow(t, 2) * 12 +
            coeffs.segment(5*n_dim, n_dim) * pow(t, 3) * 20;
      }

      return poly_interp;

    }

    static Eigen::MatrixXd quinticPolynomialEmidE(int N, // number of samples
                                                  double tf, // final time
                                                  const Eigen::VectorXd& x0, // initial state
                                                  const Eigen::VectorXd& xm, // middle state
                                                  const Eigen::VectorXd& xf, // final state
                                                  const Eigen::VectorXd& dx0, // initial velocity
                                                  const Eigen::VectorXd& dxm, // final velocity
                                                  const Eigen::VectorXd& dxf) { // final velocity


      int n_dim = x0.size();

      Eigen::MatrixXd poly_interp = Eigen::MatrixXd::Zero(3*n_dim, N);

      Eigen::VectorXd boundary_conditions_stack(6*n_dim);
      boundary_conditions_stack << x0, dx0, xm, dxm, xf, dxf;

      Eigen::MatrixXd constr_matrix = Eigen::MatrixXd::Zero(6*n_dim, 6*n_dim);

      Eigen::MatrixXd I_n = Eigen::MatrixXd::Identity(n_dim, n_dim);

      constr_matrix.block(0,0,2*n_dim,2*n_dim) = Eigen::MatrixXd::Identity(2*n_dim, 2*n_dim);

      constr_matrix.block(2*n_dim, 0, n_dim, n_dim) = I_n;
      constr_matrix.block(2*n_dim, n_dim, n_dim, n_dim) = tf/2 * I_n;
      constr_matrix.block(2*n_dim, 2*n_dim, n_dim, n_dim) = pow(tf/2, 2) * I_n;
      constr_matrix.block(2*n_dim, 3*n_dim, n_dim, n_dim) = pow(tf/2, 3) * I_n;
      constr_matrix.block(2*n_dim, 4*n_dim, n_dim, n_dim) = pow(tf/2, 4) * I_n;
      constr_matrix.block(2*n_dim, 5*n_dim, n_dim, n_dim) = pow(tf/2, 5) * I_n;

      constr_matrix.block(3*n_dim, n_dim, n_dim, n_dim) = I_n;
      constr_matrix.block(3*n_dim, 2*n_dim, n_dim, n_dim) = tf * I_n;
      constr_matrix.block(3*n_dim, 3*n_dim, n_dim, n_dim) = 3 * pow(tf/2, 2) * I_n;
      constr_matrix.block(3*n_dim, 4*n_dim, n_dim, n_dim) = 4 * pow(tf/2, 3) * I_n;
      constr_matrix.block(3*n_dim, 5*n_dim, n_dim, n_dim) = 5 * pow(tf/2, 4) * I_n;

      constr_matrix.block(4*n_dim, 0, n_dim, n_dim) = I_n;
      constr_matrix.block(4*n_dim, n_dim, n_dim, n_dim) = tf * I_n;
      constr_matrix.block(4*n_dim, 2*n_dim, n_dim, n_dim) = pow(tf, 2) * I_n;
      constr_matrix.block(4*n_dim, 3*n_dim, n_dim, n_dim) = pow(tf, 3) * I_n;
      constr_matrix.block(4*n_dim, 4*n_dim, n_dim, n_dim) = pow(tf, 4) * I_n;
      constr_matrix.block(4*n_dim, 5*n_dim, n_dim, n_dim) = pow(tf, 5) * I_n;

      constr_matrix.block(5*n_dim, n_dim, n_dim, n_dim) = I_n;
      constr_matrix.block(5*n_dim, 2*n_dim, n_dim, n_dim) = 2 * tf * I_n;
      constr_matrix.block(5*n_dim, 3*n_dim, n_dim, n_dim) = 3 * pow(tf, 2) * I_n;
      constr_matrix.block(5*n_dim, 4*n_dim, n_dim, n_dim) = 4 * pow(tf, 3) * I_n;
      constr_matrix.block(5*n_dim, 5*n_dim, n_dim, n_dim) = 5 * pow(tf, 4) * I_n;


      Eigen::VectorXd coeffs = constr_matrix.inverse() * boundary_conditions_stack;

      for(int i=0; i<N; ++i) {
        double t = tf/(N-1) * i;
        poly_interp.col(i).head(n_dim) =
            coeffs.segment(0, n_dim) +
            coeffs.segment(n_dim, n_dim) * t +
            coeffs.segment(2*n_dim, n_dim) * pow(t, 2) +
            coeffs.segment(3*n_dim, n_dim) * pow(t, 3) +
            coeffs.segment(4*n_dim, n_dim) * pow(t, 4) +
            coeffs.segment(5*n_dim, n_dim) * pow(t, 5);

        poly_interp.col(i).segment(n_dim,n_dim) =
            coeffs.segment(n_dim, n_dim) +
            coeffs.segment(2*n_dim, n_dim) * t * 2 +
            coeffs.segment(3*n_dim, n_dim) * pow(t, 2) * 3 +
            coeffs.segment(4*n_dim, n_dim) * pow(t, 3) * 4 +
            coeffs.segment(5*n_dim, n_dim) * pow(t, 4) * 5;

        poly_interp.col(i).tail(n_dim) =
            coeffs.segment(n_dim, n_dim) +
            coeffs.segment(2*n_dim, n_dim) * 2 +
            coeffs.segment(3*n_dim, n_dim) * t * 6 +
            coeffs.segment(4*n_dim, n_dim) * pow(t, 2) * 12 +
            coeffs.segment(5*n_dim, n_dim) * pow(t, 3) * 20;
      }
      return poly_interp;
    }


} // end namespace utils















