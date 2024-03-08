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
using math_util::DegToRad;
using math_util::RadToDeg;

using vector_map::VectorMap;

DEFINE_double(num_particles, 40, "Number of particles");

namespace particle_filter {

config_reader::ConfigReader config_reader_({"config/particle_filter.lua"});

ParticleFilter::ParticleFilter() :
    prev_odom_loc_(0, 0),
    prev_odom_angle_(0),
    odom_initialized_(false),
    sum_weight(1),
    num_valid_particles(FLAGS_num_particles),
    num_updates_done(0),
    num_updates_reqd_for_resample(3),
    dist_traveled(0),
    ith_ray(10),
    deg_offset(ith_ray * laser_ang_deg_res),
    debug_print(false) {}

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
                                            std::vector<Eigen::Vector2f>* scan_ptr) {

    std::vector<Eigen::Vector2f>& scan = *scan_ptr;

    Eigen::Vector2f laser_loc(0, 0);
    laser_loc.x() = loc.x() + sin(angle) * 0.2;
    laser_loc.y() = loc.y() + cos(angle) * 0.2;

    scan.resize(num_ranges / ith_ray);

    vector<float> scan_min_dists;
    for(int i = 0; i < num_ranges / ith_ray; i++) scan_min_dists.push_back(std::numeric_limits<float>::max());
        
    float degtoradresult = math_util::DegToRad(deg_offset);
    // iterate through scan to set points for each ray
    for (size_t j = 0; j < map_.lines.size(); ++j) {
        // iterate through map to check for collisions
        float alpha = angle + angle_min - degtoradresult;
        Eigen::Vector2f rm_pt(0, 0);
        for (size_t i = 0; i < scan.size(); ++i) {
            // to check for collisions, construct a line2f from range_min to range_max, in the direction of the ray, centered around laser pose
            const line2f map_line = map_.lines[j];
            // need to use math to calculate endpoint of the line based on range_max. treat laser location as point 0
            // 10 is a magic number rn describing the angle increment. we should tune that (and by extension num_ranges)
            alpha += degtoradresult;
            // the range max point
            // Eigen::Vector2f rm_pt(range_max * sin(alpha) + laser_loc.x(), );
            float cosalpha = cos(alpha);
            float sinalpha = sin(alpha);
            Eigen::Vector2f min_pt(range_min * sinalpha + laser_loc.x(), range_min * cosalpha + laser_loc.y()); 
            rm_pt.x() = range_max * cosalpha + laser_loc.x();
            rm_pt.y() = range_max * sinalpha + laser_loc.y();
            line2f my_line(min_pt.x(), min_pt.y(), rm_pt.x(), rm_pt.y());

            // check for intersection with this line and the map line
            Vector2f intersection_point; // Return variable
            bool intersects = map_line.Intersection(my_line, &intersection_point);
            
            if (intersects) {
                // if intersection exists, "first" collision wins
                float curr_dist = pow(intersection_point.x() - laser_loc.x(), 2) + pow(intersection_point.y() - laser_loc.y(), 2);
                if (curr_dist < scan_min_dists[i]) {
                    scan_min_dists[i] = curr_dist;
                    scan[i] = intersection_point;
                }
            } else if(scan_min_dists[i] == std::numeric_limits<float>::max()) {
                scan[i] = rm_pt;
                scan_min_dists[i] = std::numeric_limits<float>::max() - 1;
            }
        }
    }
}








// float calc_variance(const vector<float>& ranges) {
//     float avg = 0;
//     for(size_t i = 0; i < ranges.size(); ++i) {
//         avg += ranges[i];
//     }
//     avg /= ranges.size();
//     float variance = 0;
//     for(size_t i = 0; i < ranges.size(); ++i) {
//         variance += (ranges[i] - avg) * (ranges[i] - avg);
//     }
//     return variance /= ranges.size();
// }







void ParticleFilter::Update(const vector<float>& ranges,
                            float range_min,
                            float range_max,
                            float angle_min,
                            float angle_max,
                            Particle* p_ptr) {
    // printf("update\n");
    int num_ranges = ranges.size() / ith_ray;
    std::vector<Eigen::Vector2f> scan_ptr(num_ranges);
    for(int i = 0; i < num_ranges; i++) {
        scan_ptr[i] = Vector2f(0, 0);
    }
    
    GetPredictedPointCloud(p_ptr->loc, p_ptr->angle, ranges.size(), range_min, range_max, angle_min, angle_max, &scan_ptr);
    double log_lik = 0.0;
    float divisor = 0.05;
    // robustification will be on 14 - Expecting The Unexpected slide 28
    float d_short = 1;
    float d_short_squared = pow(d_short, 2);
    float d_long = 8;
    float d_long_squared = pow(d_long, 2);
    float delta = 0;

    for (size_t i = 0; i < scan_ptr.size(); ++i) {
        // s_hat is (dist btwn laser and scan[i] points) - range_min
        // TUNABLE: check if this should be sqnorm instead of norm if particle filter is slow
        // TODO: check if I need to subtract laser_loc.x and .y from here (and also figure out how to calc that)
        
        Vector2f s_hat_pts(scan_ptr[i] - p_ptr->loc);
        double s_hat_dist = s_hat_pts.norm();
        float s_hat = s_hat_dist - range_min;
        float s = ranges[i * ith_ray] - range_min;
        // range_min 0.020000, range_max 10.000000
        if (s < range_min || s > range_max) {
          continue;
        } else if (s < s_hat - d_short) {
          delta = d_short_squared / divisor;
        } else if (s > s_hat + d_long) {
          delta = d_long_squared / divisor;
        } else {
          delta = pow((s - s_hat), 2) / divisor;
        }

        // if(delta > 1000000) {
        //     printf("SCARILY LARGE S DELTA DETECTED?\n");
        //     printf("scan_ptr x: %f, scan_ptr y: %f, pptr x: %f, pptr y: %f, s_hat: %f, s: %f, delta: %f\n", scan_ptr[i].x(), scan_ptr[i].y(), p_ptr->loc.x(), p_ptr->loc.y(), s_hat, s, delta);
        // }

        log_lik += delta;
    }

    p_ptr->weight = log_lik * 2 / num_ranges * -1;

    // if(p_ptr->weight > 1000000) {
    //     printf("RESULTING IN CONCERNINGLY LARGE PARTICLE WEIGHT?\n");
    //     printf("weight: %f, particle x: %f, particle y: %f, particle angle: %f\n", p_ptr->weight, p_ptr->loc.x(), p_ptr->loc.y(), p_ptr->angle);
    //     printf("log_lik: %f", log_lik);
    // }

}







void ParticleFilter::Resample() {
    // printf("resample\n");
    vector<Particle> new_particles;
    if(sum_weight == 0) {
        printf("sum is zero, returning early from resample\n");
        return; 
    } 
    double step_size = sum_weight / FLAGS_num_particles;
    double new_weight = 1 / FLAGS_num_particles;
    sum_weight = 1;
    double last = rng_.UniformRandom(0, step_size) - step_size;
    int index = 0;
    double sum = 0;

    // printf("resample: sum_weight: %f, step_size: %f, new_weight: %f, last: %f, index: %d, sum: %f\n", sum_weight, step_size, new_weight, last, index, sum);

    while(new_particles.size() < FLAGS_num_particles) {
        // printf("resample outer loop\n");
        // if(particles_[index].weight == 0) {
        //     index++;
        //     continue;
        // }
        sum += particles_[index].weight;
        while(last + step_size < sum) {
            // printf("resample inner loop\n");
            Particle p;
            p.loc.x() = particles_[index].loc.x();
            p.loc.y() = particles_[index].loc.y();
            p.angle = particles_[index].angle;
            p.weight = new_weight;
            new_particles.push_back(p);
            last += step_size;
        }
        index++;
    }

    // printf("done resampling\n");
    // After resampling:
    particles_ = new_particles;
    num_valid_particles = FLAGS_num_particles;
    num_updates_done = 0;
}








void ParticleFilter::ObserveLaser(const vector<float>& ranges,
                                  float range_min,
                                  float range_max,
                                  float angle_min,
                                  float angle_max) {
    // printf("observelaser\n");
    sum_weight = 0;
    double max_likelihood = 0;
    // aka if dist_traveled < 0.387298335 meters
    // TUNABLE: the .15 thing
    if (dist_traveled < .15) return;
    dist_traveled = 0;

    for (size_t i = 0; i < FLAGS_num_particles; ++i) {
        // if (particles_[i].weight == 0) continue;
        ParticleFilter::Update(ranges, range_min, range_max, angle_min, angle_max, &(particles_[i]));
        double likelihood = particles_[i].weight;
        if(max_likelihood > likelihood) max_likelihood = likelihood;
    }
    // normalize all weights (see 16 - Problems with Particle Filters slide 10)
    for (size_t i = 0; i < FLAGS_num_particles; ++i) {
        // if(particles_[i].weight == 0) continue;
        particles_[i].weight = abs((particles_[i].weight - max_likelihood));
        // if(particles_[i].weight == 0) num_valid_particles--;
        sum_weight += particles_[i].weight;
    }
    num_updates_done++;

    if (num_updates_done >= num_updates_reqd_for_resample) ParticleFilter::Resample();
}








double magnitude(float x, float y) {
  return sqrt(pow(x, 2) + pow(y, 2));
}







void ParticleFilter::Predict(const Vector2f& odom_loc,
                             const float odom_angle) {

    // printf("[predict]\n");
    if (!odom_initialized_) {
        prev_odom_loc_.x() = odom_loc.x();
        prev_odom_loc_.y() = odom_loc.y();
        prev_odom_angle_ = odom_angle;
        odom_initialized_ = true;
        return;
    }

    double x_hat = odom_loc.x() - prev_odom_loc_.x();
    double y_hat = odom_loc.y() - prev_odom_loc_.y();
    double theta_hat = odom_angle - prev_odom_angle_;

    Eigen::Rotation2Df r_prev_odom(-prev_odom_angle_);

    if(x_hat == 0 && y_hat == 0 && theta_hat == 0) {
        // printf("no movement, returning early\n");
        return;
    }

    Vector2f T_odom(x_hat, y_hat);
    // squared version
    dist_traveled += pow(x_hat, 2) + pow(y_hat, 2);

    Vector2f T_delta_bl = r_prev_odom * T_odom;
    // printf("[predict] t_delta_bl x: %f, y: %f, norm: %f\n", T_delta_bl.x(), T_delta_bl.y(), T_delta_bl.norm());

    // printf("\n\n\n\nT_delta_bl: x: %f, y: %f, theta: %f\n", T_delta_bl.x(), T_delta_bl.y(), odom_angle);
    
    // tunable parameters
    float k_1 = 0.5;  //   x,y stddev's   mag(x,y) weight
    float k_2 = 0.05;  //   x,y stddev's   mag(theta) weight
    float k_3 = 0.1;   // theta stddev's   mag(x,y) weight
    float k_4 = 0.05; // theta stddev's   mag(theta) weight

    // float particle_x_hat = 0;
    // float particle_y_hat = 0;
    // float particle_theta_hat = 0;

    for (size_t i = 0; i < FLAGS_num_particles; ++i){
        // if(particles_[i].weight == 0) continue;
        Eigen::Rotation2Df r_map(particles_[i].angle);
        Vector2f T_map_one(particles_[i].loc.x(), particles_[i].loc.y());

        double xy_stddev = k_1*T_delta_bl.norm() + k_2*abs(theta_hat);
        double theta_stddev = k_3*T_delta_bl.norm() + k_4*abs(theta_hat);

        float e_x = 0;
        float e_y = 0;
        float e_theta = 0;

        if(xy_stddev != 0) {
            e_x = rng_.Gaussian(0.0, xy_stddev);
            e_y = rng_.Gaussian(0.0, xy_stddev);
        }
        if(theta_stddev != 0) {
            e_theta = rng_.Gaussian(0.0, theta_stddev);
        }

        Vector2f noise(e_x, e_y);
        Vector2f T_map = T_map_one + r_map * (T_delta_bl + noise);

        // printf("particle[%ld] old x: %f -> new x: %f\n", i, particles_[i].loc.x(), T_map.x());
        // printf("particle[%ld] old y: %f -> new y: %f\n", i, particles_[i].loc.y() , T_map.y());
        // printf("particle[%ld] old theta: %f -> new theta: %f\n\n", i, particles_[i].angle, particles_[i].angle + theta_hat + e_theta);

        particles_[i].loc.x() = T_map.x();
        particles_[i].loc.y() = T_map.y();
        particles_[i].angle += theta_hat + e_theta;

        // particle_x_hat += (r_map * T_delta_bl).x() + e_x;
        // particle_y_hat += (r_map * T_delta_bl).y() + e_y;
        // particle_theta_hat += theta_hat + e_theta;
        
        // line2f particle_line(T_map_one.x(), T_map_one.y(), particles_[i].loc.x(), particles_[i].loc.y());

        // for (size_t j = 0; j < map_.lines.size(); ++j) {
        //     const line2f map_line = map_.lines[j];
        //     Vector2f intersection_point;
        //     bool intersects = map_line.Intersection(particle_line, &intersection_point);
        //     if (intersects) {
        //         particles_[i].weight = 0;
        //         num_valid_particles--;
        //     } 
        // }
    }

    // printf("\n particle x hat: %f, actual odom x hat: %f\n");
    // printf("particle x hat: %f, actual odom x hat: %f\n");
    // printf("particle x hat: %f, actual odom x hat: %f\n");

    prev_odom_loc_ = odom_loc;
    prev_odom_angle_ = odom_angle;
}







void ParticleFilter::Initialize(const string& map_file,
                                const Vector2f& loc,
                                const float angle) {
    // printf("INITIALIZED!\n");

    map_.Load(map_file);
    vector<Particle> new_particles;

    if(debug_print) printf("INITIALIZE DATA: x: %f, y: %f, angle: %f\n", loc.x(), loc.y(), angle);

    num_valid_particles = FLAGS_num_particles;
    odom_initialized_ = false;
    num_updates_done = 0;
    dist_traveled = 0;
    double new_weights = 1 / FLAGS_num_particles;
    // initialize vector of particles with GetParticles
    for (size_t i = 0; i < FLAGS_num_particles; ++i){
        Particle p = Particle();
        // init in Gaussian distribution around loc and angle
        p.loc.x() = loc.x() + rng_.Gaussian(0.0, 0.001);
        p.loc.y() = loc.y() + rng_.Gaussian(0.0, 0.001);
        p.angle = angle + rng_.Gaussian(0.0, 0.02);
        p.weight = new_weights;
        sum_weight += p.weight;
        new_particles.push_back(p);
    }

    particles_ = new_particles;
}







void ParticleFilter::GetLocation(Eigen::Vector2f* loc_ptr, 
                                 float* angle_ptr) const {
    // printf("getlocation\n");
    Vector2f& loc = *loc_ptr;
    float& angle = *angle_ptr;

    double x_locs = 0.0;
    double y_locs = 0.0;
    double sines = 0.0;
    double cosines = 0.0;

    double local_sum_weight = 0;

    for (size_t i = 0; i < FLAGS_num_particles; ++i){
        if(particles_[i].weight == 0) continue;
        // printf("getlocation -- particle[%ld] -- weight: %f\n", i, particles_[i].weight);
        x_locs += particles_[i].loc.x() * particles_[i].weight;
        y_locs += particles_[i].loc.y() * particles_[i].weight;
        sines += sin(particles_[i].angle) * particles_[i].weight;
        cosines += cos(particles_[i].angle) * particles_[i].weight;
        local_sum_weight += particles_[i].weight;
    }

    if(sum_weight == 0) {
        printf("this location is probably wrong, sum_weight==0\n");
        loc.x() = x_locs / 1;
        loc.y() = y_locs / 1;
        angle = atan2(sines / 1, cosines / 1);
    } else {
        loc.x() = x_locs / local_sum_weight;
        loc.y() = y_locs / local_sum_weight;
        angle = atan2(sines / local_sum_weight, cosines / local_sum_weight) ;
    }
}


}  // namespace particle_filter
