//
// Created by ataka on 2.03.2026.
//

#ifndef SAT_PERFORM_ORBITCONFIG_H
#define SAT_PERFORM_ORBITCONFIG_H

#include <cmath>

namespace config {
    struct OrbitConfig {
        double altitude; //km
        double inclination;
        double eccentricity;
        double RAAN_deg;
        double arpg_deg;
        double nu_deg;

        double semi_major_axis;
        double mean_motion;
        double orbital_period;

        static OrbitConfig create(double altitude_km, double inclination_deg) {
            OrbitConfig orb;
            orb.altitude = altitude_km;
            orb.inclination = inclination_deg;
            orb.eccentricity = 0.01;
            orb.RAAN_deg = 0.01;
            orb.nu_deg = 0.0;
            orb.arpg_deg = 270;

            constexpr double R_earth = 6378137.0;
            constexpr double mu =  3.986004418e14;
            orb.semi_major_axis = R_earth + altitude_km * 1000.0;
            orb.mean_motion = std::sqrt(mu / std::pow(orb.semi_major_axis,3));
            orb.orbital_period = 2.0 * M_PI / orb.mean_motion;

            return orb;

        }
    };
}

#endif //SAT_PERFORM_ORBITCONFIG_H