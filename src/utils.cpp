#include <utils.hpp>

// angleSignedDistance using rotation matrices
/*double utils::angleSignedDistance(double a, double b){

  Eigen::Rotation2D Ra = Eigen::Rotation2D(a);
  Eigen::Rotation2D Rb = Eigen::Rotation2D(b);

  Eigen::Rotation2D R = Rb.inverse()*Ra;
  return R.smallestAngle();
}*/

double utils::angleSignedDistance(double a, double b){
  double d = a - b;
  while(d >  M_PI) d = d - 2.0*M_PI;
  while(d < -M_PI) d = d + 2.0*M_PI;

  return d;
}

// Angular errors computed with RPY
/*Eigen::Vector3d utils::angleSignedDistance(Eigen::Vector3d a, Eigen::Vector3d b) {
  Eigen::Vector3d diff;
  diff(0) = utils::angleSignedDistance(a(0),b(0));
  diff(1) = utils::angleSignedDistance(a(1),b(1));
  diff(2) = utils::angleSignedDistance(a(2),b(2));
  return diff;
}*/

// Angular errors computed with AxisAngle
Eigen::Vector3d utils::angleSignedDistance(Eigen::Vector3d a, Eigen::Vector3d b) {
  Eigen::Matrix3d Ra = utils::rot(a);
  Eigen::Matrix3d Rb = utils::rot(b);

  Eigen::Matrix3d Rdiff = Rb.transpose() * Ra;
  auto aa = Eigen::AngleAxisd(Rdiff);
  return aa.angle() * Ra * aa.axis();
}

Eigen::Vector3d utils::getRPY(const Eigen::Matrix3d & rotMatrix) {
  Eigen::Vector3d RPY;
  RPY << atan2(rotMatrix(2, 1), rotMatrix(2, 2)),
      atan2(-rotMatrix(2, 0), sqrt(rotMatrix(2, 1) * rotMatrix(2, 1) + rotMatrix(2, 2) * rotMatrix(2, 2))),
      atan2(rotMatrix(1, 0), rotMatrix(0, 0));
  return RPY;
}

Eigen::Vector3d utils::getRPY(dart::dynamics::BodyNode *body) {
  Eigen::Matrix3d rotMatrix = body->getTransform().rotation();
  return utils::getRPY(rotMatrix);
}

Eigen::Matrix3d utils::rot(const Eigen::Vector3d & eul){

  Eigen::Matrix3d r;

  Eigen::Matrix3d rx = utils::rotx(eul(0));
  Eigen::Matrix3d ry = utils::roty(eul(1));
  Eigen::Matrix3d rz = utils::rotz(eul(2));

  r = rz*ry*rx;

  return r;
}

Eigen::Matrix2d utils::rot(const double angle){

  Eigen::Matrix2d R;
  const double c = cos(angle);
  const double s = sin(angle);
  R(0,0) = c;
  R(0,1) = -s;
  R(1,0) = s;
  R(1,1) = c;

  return R;
}



Eigen::Matrix4d utils::v2t(const Eigen::VectorXd & v){

  Eigen::Matrix4d m = Eigen::Matrix4d::Identity();

  Eigen::Vector3d eul = v.head(3);
  Eigen::Matrix3d r = utils::rot(eul);

  m.block<3,3>(0,0) = r;
  m.block<3,1>(0,3) = v.tail(3);

  return m;
}

Eigen::VectorXd utils::t2v(const Eigen::Matrix4d & m) {
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
inline Eigen::VectorXd utils::vvRel(const Eigen::VectorXd & v2, const Eigen::VectorXd & v1) {
  return utils::t2v(utils::v2t(v1).inverse()*utils::v2t(v2));
}
