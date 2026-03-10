//
// Created by ataka on 2.03.2026.
//

#ifndef SAT_PERFORM_SATELLITECONFIG_H
#define SAT_PERFORM_SATELLITECONFIG_H


#include <Eigen/Dense>

namespace config {
    enum class SatelliteType {
        Cubesat3U,
        UPMSAT2
    };

    struct SatelliteConfig {
        double mass;
        Eigen::Vector3d dimensions;
        Eigen::Matrix3d inertia;
        double max_dipole;
        bool use_dualspin;
        double I_rotor;
        double omega_rotor;

        static SatelliteConfig create(SatelliteType type, bool dual_spin) {
            SatelliteConfig cfg;
            cfg.use_dualspin = dual_spin;

            if (type == SatelliteType::Cubesat3U) {
                cfg.mass = 4;
                cfg.dimensions = {0.1,0.1,0.3};
                cfg.max_dipole = 0.2;
            } else {
                cfg.mass = 50;
                cfg.dimensions = {0.5,0.5,0.6};
                cfg.max_dipole = 15;
            }


            //Calculate inertia
            double l = cfg.dimensions(0);
            double w = cfg.dimensions(1);
            double h = cfg.dimensions(2);
            double m = cfg.mass ;
            double Ix = m/12.0 * (w*w + h*h);
            double Iy = m/12.0 * (l*l + h*h);
            double Iz = m/12.0 * (l*l + w*w);
            cfg.inertia = Eigen::Matrix3d::Identity();
            cfg.inertia(0,0) = Ix;
            cfg.inertia(1,1) = Iy;
            cfg.inertia(2,2) = Iz;

            cfg.I_rotor = dual_spin ? 0.005 : 0.0;
            cfg.omega_rotor = dual_spin ? 2*M_PI : 0.0;
            return cfg;
        }
    };
}


#endif //SAT_PERFORM_SATELLITECONFIG_H