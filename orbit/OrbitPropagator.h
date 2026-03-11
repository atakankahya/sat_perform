//
// Created by ataka on 3.03.2026.
//

#ifndef SAT_PERFORM_ORBITPROPAGATOR_H
#define SAT_PERFORM_ORBITPROPAGATOR_H

#include "Eigen/Dense"
#include "../config/OrbitConfig.h"

namespace orbit {
    // ECI Position and velocity
    struct OrbitState {
        Eigen::Vector3d r_eci;
        Eigen::Vector3d v_eci;
    };

    struct LLA {
        double lat_deg;
        double lon_deg;
        double alt_km;
    };


    class OrbitPropagator {
        /**
            * @param cfg       Orbital configuration (semi-major axis, eccentricity, etc.)
            * @param gmst0_rad GMST at simulation start [rad]. Default 0 means Greenwich
            *                  is aligned with the vernal equinox at t = 0.
            */

    public:
        explicit OrbitPropagator(const config::OrbitConfig& cfg, double gmst0_rad = 0.0);

        void propagate(double dt); // advance the orbit by dt


        [[nodiscard]] OrbitState getECIState() const;
        [[nodiscard]] LLA getLLA() const;
        [[nodiscard]] double getTrueAnomaly_deg() const;
        [[nodiscard]] double getElapsedTime() const{return elapsed_time_;};

        [[nodiscard]] double computeGMST() const;

    private:
        double a_;
        double e_;
        double inc_;
        double RAAN_;
        double argp_;
        double n_;

        double M_; //current mean anomaly
        double elapsed_time_;
        double gmst0_;

        static constexpr double mu_          = 3.986004418e14;    // Earth grav. param [m^3/s^2]
        static constexpr double R_earth_     = 6378137.0;         // WGS-84 semi-major axis [m]
        static constexpr double f_           = 1.0 / 298.257223563; // WGS-84 flattening
        static constexpr double omega_earth_ = 7.2921150e-5;      // Earth rotation rate [rad/s]

        [[nodiscard]] double solveKepler(double M, double e) const;
        [[nodiscard]] double eccentricToTrue(double E) const;
        [[nodiscard]] double trueToMean(double nu) const;
        [[nodiscard]] Eigen::Matrix3d perifocalToECI_DCM() const;
        [[nodiscard]] Eigen::Vector3d eciToECEF(const Eigen::Vector3d& r_eci) const;
        [[nodiscard]] LLA ecefToLLA(const Eigen::Vector3d& r_ecef) const;

    };
}


#endif //SAT_PERFORM_ORBITPROPAGATOR_H