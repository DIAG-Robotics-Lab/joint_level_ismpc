#pragma once

#include <Eigen/Geometry>
#include <Eigen/Core>

#include <iostream>
#include <fstream>
#include <vector>

#include <dart/dart.hpp>


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


namespace utils {

// Contruct 2D and 3D rotation matrix
Eigen::Matrix2d rot(const double angle);
Eigen::Matrix3d rot(const Eigen::Vector3d & eul);

double angleSignedDistance(double a, double b);
Eigen::Vector3d angleSignedDistance(Eigen::Vector3d a, Eigen::Vector3d b);
//Eigen::VectorXd angleSignedDistance(Eigen::VectorXd a, Eigen::VectorXd b);

// Check if an element is in a vector
template<typename T>
bool is_in_vector(const std::vector<T> & vector, const T & elt) {
  return vector.end() != std::find(vector.begin(),vector.end(),elt);
}

// Transform a BodyNode or a rotation matrix in a vector of RPY angles
Eigen::Vector3d getRPY(dart::dynamics::BodyNode *body);
Eigen::Vector3d getRPY(const Eigen::Matrix3d & rotMatrix);

Eigen::Matrix4d v2t(const Eigen::VectorXd & v);
Eigen::VectorXd t2v(const Eigen::Matrix4d & m);
// Express v2 in the frame of v1
Eigen::VectorXd vvRel(const Eigen::VectorXd & v2, const Eigen::VectorXd & v1);

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


/*class ContainerBase{
 public:
  virtual void log() = 0;
  //virtual ~ContainerBase();
};

template <typename Derived>
class Container : public ContainerBase {
  public:
  Container(Derived * addr_to_be_logged, std::string name) {
    var_ptr = addr_to_be_logged;
    createFile(name);
  }
  void log() {
    file << var_ptr->transpose() << std::endl;
  }
 private:
  Derived * var_ptr;
  std::ofstream file;

  void createFile(std::string name) {
    std::cout << "created LOG in " << realpath("../data/", NULL) + ("/" + name) << std::endl;
    file = std::ofstream(realpath("../data/", NULL) + ("/" + name), std::ofstream::out);
  }
};

class Logger {
 public:
  void add(utils::ContainerBase* var_to_be_logged) {
    list.push_back(var_to_be_logged);
  }
  void logAll() {
  for(auto logger: list) (*logger).log();
  }
 private:
  std::vector<utils::ContainerBase*> list;
};*/

} // end namespace utils















