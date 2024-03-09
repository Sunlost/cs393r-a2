#pragma once

#include "navigation_params.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include "shared/math/math_util.h"

namespace motion_primitives {
float run1DTOC(const navigation::MotionLimits &limits, const float x0, const float v0,
               const float xf, const float vf,
               const float dt = 0.05);
}