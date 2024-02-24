//========================================================================
//  This software is free: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License Version 3,
//  as published by the Free Software Foundation.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  Version 3 in the file COPYING that came with this distribution.
//  If not, see <http://www.gnu.org/licenses/>.
//========================================================================
/*!
\file    particle-filter.cc
\brief   Particle Filter Starter Code
\author  Joydeep Biswas, (C) 2019
*/
//========================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include "gflags/gflags.h"
#include "glog/logging.h"
#include "shared/math/geometry.h"
#include "shared/math/line2d.h"
#include "shared/math/math_util.h"
#include "shared/util/timer.h"

#include "config_reader/config_reader.h"
#include "particle_filter.h"

#include "vector_map/vector_map.h"

using geometry::line2f;
using std::cout;
using std::endl;
using std::string;
using std::swap;
using std::vector;
using Eigen::Vector2f;
using Eigen::Vector2i;
using vector_map::VectorMap;

DEFINE_double(num_particles, 50, "Number of particles");

namespace particle_filter {

config_reader::ConfigReader config_reader_({"config/particle_filter.lua"});

ParticleFilter::ParticleFilter() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false) {}

void ParticleFilter::GetParticles(vector<Particle>* particles) const {
  *particles = particles_;
}

void ParticleFilter::GetPredictedPointCloud(const Vector2f& loc,
                                            const float angle,
                                            int num_ranges,
                                            float range_min,
                                            float range_max,
                                            float angle_min,
                                            float angle_max,
                                            vector<Vector2f>* scan_ptr) {
  // loc // robot's pose loc
  // angle // robot's pose angle
  // num_ranges // number of rays to use
  // range_min // Minimum observable range
  // range_max // Maximum observable range
  // angle_min // Angle of the last ray
  // angle_max // Angle of the first ray
  // scan_ptr // save the predicted point cloud? aka save point of either intersection with wall or point at range_max

  vector<Vector2f>& scan = *scan_ptr;
  // Compute what the predicted point cloud would be, if the car was at the pose
  // loc, angle, with the sensor characteristics defined by the provided
  // parameters.
  // This is NOT the motion model predict step: it is the prediction of the
  // expected observations, to be used for the update step.

  // Note: The returned values must be set using the `scan` variable:
  scan.resize(num_ranges);
  // Fill in the entries of scan using array writes, e.g. scan[i] = ...
  for (size_t i = 0; i < scan.size(); ++i) {
    scan[i] = Vector2f(0, 0);
  }

  // The line segments in the map are stored in the `map_.lines` variable. You
  // can iterate through them as:
  for (size_t i = 0; i < map_.lines.size(); ++i) {
    const line2f map_line = map_.lines[i];
    // The line2f class has helper functions that will be useful.
    // You can create a new line segment instance as follows, for :
    line2f my_line(1, 2, 3, 4); // Line segment from (1,2) to (3.4).
    // Access the end points using `.p0` and `.p1` members:
    printf("P0: %f, %f P1: %f,%f\n", 
           my_line.p0.x(),
           my_line.p0.y(),
           my_line.p1.x(),
           my_line.p1.y());

    // Check for intersections:
    bool intersects = map_line.Intersects(my_line);
    // You can also simultaneously check for intersection, and return the point
    // of intersection:
    Vector2f intersection_point; // Return variable
    intersects = map_line.Intersection(my_line, &intersection_point);
    if (intersects) {
      printf("Intersects at %f,%f\n", 
             intersection_point.x(),
             intersection_point.y());
    } else {
      printf("No intersection\n");
    }
  }

  // TODO:
  // iterate through scan to set points as in line 91
  // iterate through map to check for collisions as in line 97-114
    // to check for collisions, construct a line2f from range_min to range_max, in the direction of the ray, centered around laser pose
    // check for intersection with this line and the map line
    // if intersection exists, "first" collision wins, so use continue and go on to next 
    // else if no collision, set scan[i] to the point at range_max

}

void ParticleFilter::Update(const vector<float>& ranges,
                            float range_min,
                            float range_max,
                            float angle_min,
                            float angle_max,
                            Particle* p_ptr) {

  // Implement the update step of the particle filter here.
  // You will have to use the `GetPredictedPointCloud` to predict the expected
  // observations for each particle, and assign weights to the particles based
  // on the observation likelihood computed by relating the observation to the
  // predicted point cloud.

  // TODO: STEP 1
  // use GetPredictedPointCloud to predict expected observations for particle conditioned on the map
    // init a scan_ptr vector of points... size should be num_ranges
    // pass in robot's location and angle for loc and angle
    // num_ranges should equal something like (angle_max - angle_min) / 10, so every 10 degrees we use the lidar range
    // ParticleFilter::GetPredictedPointCloud(const Vector2f& loc,
    //                                           const float angle,
    //                                           int num_ranges,
    //                                           range_min,
    //                                           range_max,
    //                                           angle_min,
    //                                           angle_max,
    //                                           vector<Vector2f>* scan_ptr)

  // TODO: STEP 2
  // compare particle observation to prediction
  // ranges holds our observation s, scan_ptr (kinda) holds our prediction s_hat
  // loop over both of these
    // tunable param: sd
    // s_hat is basically (dist btwn laser and scan[i] points) - range_min
    // but what is s
    // anyway
    // log_lik -= ((s_hat - s)/sd^2)^2

  // TODO: STEP 3
  // assign weight to particle
  // particle struct has field weight, set that to log_lik?
  /* Dr. Biswas's answer from 2/14
     These [point-cloud point location] estimations are based on the particles that we are modeling, 
     trying to guess what the point cloud would look like from that particle's perspective 
     
     If the particle's predicted point_cloud matches with the "ground truth" (what is this, ask Amanda)
     that is considered a "good" particle, and "bad" vice-versa. This will inform our weightings on the 
     particles. 
     
     need more clarification on what this "ground truth" is and how we use it...
  */
}

void ParticleFilter::Resample() {
  // Resample the particles, proportional to their weights.
  // The current particles are in the `particles_` variable. 
  // Create a variable to store the new particles, and when done, replace the
  // old set of particles:
  // vector<Particle> new_particles';
  // During resampling: 
  //    new_particles.push_back(...)
  // After resampling:
  // particles_ = new_particles;

  /*
    ED: Either we keep some tunable absolute threshold for the weights for discarding
    old unlikely particles 
    OR
    could figure out the range of of weights we have just calculated, figure out threshold and remove particles. 
    Don't filter out too many particles, maybe start by dropping 30% of the range of weights? 
    
    The latter sounds more correct.

    also need to clarify what they mean about "duplication", whether that is copying over an old point 
    from the old particles to the new ones 
    or 
    actually duplicating the particles s.t where there was one particle present in the old set, now there would be two. 
    
    The former sounds like the correct way.
  */ 

  // You will need to use the uniform random number generator provided. For
  // example, to generate a random number between 0 and 1:
  float x = rng_.UniformRandom(0, 1);
  printf("Random number drawn from uniform distribution between 0 and 1: %f\n",
         x);
}

void ParticleFilter::ObserveLaser(const vector<float>& ranges,
                                  float range_min,
                                  float range_max,
                                  float angle_min,
                                  float angle_max) {
  // A new laser scan observation is available (in the laser frame)
  // Call the Update and Resample steps as necessary.

  // TODO STEP 1: figure out how to call update
  // init a num_updates variable
  // tunable param: d
  // if we have traveled at least distance d
      // loop over particle vector
        // ParticleFilter::Update(ranges, range_min, range_max, angle_min, angle_max,
        //                           Particle* p_ptr)
        // increment num_updates

  // TODO STEP 2: figure out how to call resample
  // tunable param: n
  // if it has been n updates since our last resample
    // ParticleFilter::Resample()
}

void ParticleFilter::Predict(const Vector2f& odom_loc,
                             const float odom_angle) {
  // Implement the predict step of the particle filter here.
  // A new odometry value is available (in the odom frame)
  // Implement the motion model predict step here, to propagate the particles
  // forward based on odometry.


  // You will need to use the Gaussian random number generator provided. For
  // example, to generate a random number from a Gaussian with mean 0, and
  // standard deviation 2:
  float x = rng_.Gaussian(0.0, 2.0);
  printf("Random number drawn from Gaussian distribution with 0 mean and "
         "standard deviation of 2 : %f\n", x);

  // inside of a loop over particles
  // TODO: STEP 1
  // tunable params: k_1, k_2, k_3, k_4
  // calc x_hat, y_hat, and theta_hat from odom_loc and odom_angle

  // generate e_x and e_y from a distribution with a standard deviation of k_1*mag(x_hat^2 + y_hat^2) + k_2*abs(theta_hat)

  // generate e_theta from a distribution with a standard deviation of k_3*mag(x_hat^2 + y_hat^2) + k_4*abs(theta_hat)

  // d_x = x_hat + e_x
  // d_y = y_hat + e_y
  // d_theta = theta_hat + e_theta

  // TODO: STEP 2
  // calculate robot's current location with GetLocation()
  // model a particle some d_x, d_y, and d_theta away from current location
  // push particle to particle vector

}

void ParticleFilter::Initialize(const string& map_file,
                                const Vector2f& loc,
                                const float angle) {
  // The "set_pose" button on the GUI was clicked, or an initialization message
  // was received from the log. Initialize the particles accordingly, e.g. with
  // some distribution around the provided location and angle.
  map_.Load(map_file);

  // TODO
  // initialize vector of particles with GetParticles
  // initialize location estimate for robot

}

void ParticleFilter::GetLocation(Eigen::Vector2f* loc_ptr, 
                                 float* angle_ptr) const {
  Vector2f& loc = *loc_ptr;
  float& angle = *angle_ptr;
  // Compute the best estimate of the robot's location based on the current set
  // of particles. The computed values must be set to the `loc` and `angle`
  // variables to return them. Modify the following assignments:
  loc = Vector2f(0, 0);
  angle = 0;

  // TODO
  // this is different from mean pose -- need to figure out where that goes, maybe in predict?
  // pick particle based on the best observation likelihood (aka its weight)
}


}  // namespace particle_filter
