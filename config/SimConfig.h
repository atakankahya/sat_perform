//
// Created by ataka on 2.03.2026.
//

#ifndef SAT_PERFORM_SIMCONFIG_H
#define SAT_PERFORM_SIMCONFIG_H
#include "Eigen/Dense"
#include "../utils/Conversions.h"


namespace config {
    struct SimConfig {
        //Simulation configuration
        double dt = 0.1;
        double control_period = 2.0;
        double t_episode = 8000.0;
        Eigen::Vector3d omega_init = {-0.1, 0.1, 0};
        Eigen::Vector3d omega_desired = {0, 0, 0};

        //Initial attitude
        double roll_deg = 0;
        double pitch_deg = 30;
        double yaw_deg = 60;

        //Start date
        std::array<double, 6> start_date = {2026,3,10,0,0,0};
    };
}

#endif //SAT_PERFORM_SIMCONFIG_H