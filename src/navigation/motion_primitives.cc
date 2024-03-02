#include "motion_primitives.h"

#include <iostream>

using std::cout;
using std::endl;
using std::max;
using std::min;

namespace motion_primitives {

float getStoppingDistance(const float v, const navigation::MotionLimits &limits) {
  return math_util::Sq(v) / (2 * limits.max_deceleration);  // return meters
}

float run1DTOC(const navigation::MotionLimits &limits,
               const float x0,
               const float v0,
               const float xf,
               const float vf,
               const float dt  // timestep
) {
  // Note This function only works with positive velocities
  char phase = '?';  // A: accelerate, C: cruise, D: deccelerate, X: stop

  float dv_a = limits.max_acceleration * dt;   // in m/s
  float dv_d = limits.max_deceleration * dt;  // in m/s

  const float distance_left = xf - x0;
  const float v0_stopping_distance = getStoppingDistance(v0, limits);
  if (distance_left <= 0 || v0_stopping_distance > distance_left) {
    phase = 'X';
    cout << "WARNING!! OVERSHOOT" << endl;
    cout << "distance_left: " << distance_left << ", v0_stopping_distance"
         << v0_stopping_distance << endl;

    return max<float>(0, v0 - dv_d);
  }

  float v_accel = min(v0 + dv_a, limits.max_speed);

  float cruise_stopping_distance = v0 * dt + v0_stopping_distance;
  float accelerate_stopping_distance =
      0.5 * (v0 + v_accel) * dt + getStoppingDistance(v_accel, limits);

  float toc_velocity = 0;
  if (distance_left > accelerate_stopping_distance) {
    // Check if we have enough distance if we accelerate
    phase = 'A';
    toc_velocity = v_accel;
  } else if (distance_left > cruise_stopping_distance) {
    // Check if we have enough distance if we cruise
    phase = 'C';
    toc_velocity = v0;
  } else {
    // If not enough stopping distance at current speed
    phase = 'D';
    toc_velocity = max<float>(v0 - dv_d, 0.0);
  }

  // Add debugging printout here
  // bool debug = true;
  if (true) {
    std::ofstream file("1DTOC.txt");  // Using constructor to open the file
    printf("\nPhase: %c\n", phase);
    file << "x0 " << x0 << " v0 " << v0 << " xf " << xf << " vf " << vf << " toc_v"
         << toc_velocity << '\n';
    file.close();
  }

  return toc_velocity;
}

}  // namespace motion_primitives
