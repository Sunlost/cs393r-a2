#include "constant_curvature_arc.h"

#include "motion_primitives.h"
#include "shared/math/math_util.h"

using std::fabs;
using math_util::Sign;

namespace {
const float kEpsilon = 1e-5;
}  // namespace

namespace motion_primitives {
void ConstantCurvatureArc::getControlOnCurve(
    const navigation::MotionLimits& linear_limits,
    const float linear_vel,
    const float dt,
    float& cmd_linear_vel) {

  cmd_linear_vel =
      run1DTOC(linear_limits, 0, linear_vel, arc_length_, cmd_linear_vel, dt);
}

Eigen::Vector2f ConstantCurvatureArc::getEndPoint() {
  /**
   * @brief Computes the xy endpoint of the arc in the local base_link frame
   * 
   */
  if (std::abs(curvature_) < kEpsilon) {
    return Eigen::Vector2f(arc_length_, 0);
  }
  
  const float r = 1 / curvature_;
  const float theta = arc_length_ * curvature_;
  return Eigen::Vector2f(r * std::sin(theta), r * (1 - std::cos(theta)));
}
}  // namespace motion_primitives
