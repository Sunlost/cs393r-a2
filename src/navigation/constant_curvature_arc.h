#ifndef CONSTANT_CURVATURE_ARC_H
#define CONSTANT_CURVATURE_ARC_H

#include "eigen3/Eigen/Dense"
#include "motion_primitives.h"

namespace motion_primitives {
class ConstantCurvatureArc {
 public:
  ConstantCurvatureArc() : curvature_(0), arc_length_(0), clearance_(0) {}
  ConstantCurvatureArc(float curvature)
      : curvature_(curvature), arc_length_(0), clearance_(0) {}
  ConstantCurvatureArc(float curvature, float arc_length)
      : curvature_(curvature), arc_length_(arc_length) {}
  ConstantCurvatureArc(float curvature, float arc_length, float clearance)
      : curvature_(curvature), arc_length_(arc_length), clearance_(clearance) {}

  void getControlOnCurve(const navigation::MotionLimits& linear_limits,
                         const float linear_vel,
                         const float dt,
                         float& cmd_linear_vel);

  void set_arc_length(float arc_length) { arc_length_ = arc_length; }

  void set_clearance(float clearance) { clearance_ = clearance; }

  Eigen::Vector2f getEndPoint();

  float curvature() const { return curvature_; }

  float arc_length() const { return arc_length_; }

  float clearance() const { return clearance_; }

 private:
  float curvature_;
  float arc_length_;
  float clearance_;
};
}  // namespace motion_primitives

#endif  // CONSTANT_CURVATURE_ARC_H