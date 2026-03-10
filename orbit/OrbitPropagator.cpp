//
// Created by ataka on 3.03.2026.
//

#include "OrbitPropagator.h"

#include <cmath>
#include <stdexcept>


namespace orbit {
    OrbitPropagator::OrbitPropagator(const config::OrbitConfig &cfg, double gmst0_rad)
        : a_(cfg.semi_major_axis)
        , e_(cfg.eccentricity)
        , inc_(cfg.inclination * M_PI / 180.0)
        , RAAN_(cfg.RAAN_deg * M_PI / 180.0)
        , argp_(cfg.arpg_deg * M_PI / 180.0)
        , n_(cfg.mean_motion)
        , M_(0.0)
        , elapsed_time_(0.0)
        , gmst0_(gmst0_rad)
    {
        double nu0 = cfg.nu_deg * M_PI / 180.0;
        M_ = trueToMean(nu0);
    }


    void OrbitPropagator::propagate(double dt) {
        M_ += n_ * dt;
        M_ = std::fmod(M_, 2.0 * M_PI);
        if (M_ < 0.0) {
            M_ += 2.0 * M_PI;
        }
        elapsed_time_ += dt;
    }

    double OrbitPropagator::solveKepler(double M, double e) const {
        double E = (e > 0.8) ? M_PI : M;
        for (int iter = 0; iter < 100; ++iter) {
            double f  = E - e * std::sin(E) - M;
            double fp = 1.0 - e * std::cos(E);
            double dE = f / fp;
            E -= dE;
            if (std::abs(dE) < 1e-12) break;
        }
        return E;
    }


    double OrbitPropagator::eccentricToTrue(double E) const {
        return 2.0 * std::atan2(
            std::sqrt(1.0 + e_) * std::sin(E / 2.0),
            std::sqrt(1.0 - e_) * std::cos(E / 2.0));
    }

    double OrbitPropagator::trueToMean(double nu) const {
        double E = 2.0 * std::atan2(
            std::sqrt(1.0 - e_) * std::sin(nu / 2.0),
            std::sqrt(1.0 + e_) * std::cos(nu / 2.0));
        double M = E - e_ * std::sin(E);
        if (M < 0.0) M += 2.0 * M_PI;
        return M;
    }


    double OrbitPropagator::getTrueAnomaly_deg() const {
        double E  = solveKepler(M_, e_);
        double nu = eccentricToTrue(E);
        return nu * 180.0 / M_PI;
    }

    OrbitState OrbitPropagator::getECIState() const {
        double E  = solveKepler(M_, e_);
        double nu = eccentricToTrue(E);

        double r = a_ * (1.0 - e_ * std::cos(E));

        Eigen::Vector3d r_pf(r * std::cos(nu), r * std::sin(nu), 0.0);

        double p     = a_ * (1.0 - e_ * e_);
        double coeff = std::sqrt(mu_ / p);
        Eigen::Vector3d v_pf(-coeff * std::sin(nu),
                              coeff * (e_ + std::cos(nu)),
                              0.0);

        Eigen::Matrix3d R = perifocalToECI_DCM();

        OrbitState st;
        st.r_eci = R * r_pf;
        st.v_eci = R * v_pf;
        return st;
    }

    LLA OrbitPropagator::getLLA() const {
        OrbitState st          = getECIState();
        Eigen::Vector3d r_ecef = eciToECEF(st.r_eci);
        return ecefToLLA(r_ecef);
    }


    Eigen::Matrix3d OrbitPropagator::perifocalToECI_DCM() const {
        double cO = std::cos(RAAN_), sO = std::sin(RAAN_);
        double ci = std::cos(inc_),  si = std::sin(inc_);
        double cw = std::cos(argp_), sw = std::sin(argp_);

        Eigen::Matrix3d R;
        R(0, 0) =  cO * cw - sO * sw * ci;
        R(0, 1) = -cO * sw - sO * cw * ci;
        R(0, 2) =  sO * si;

        R(1, 0) =  sO * cw + cO * sw * ci;
        R(1, 1) = -sO * sw + cO * cw * ci;
        R(1, 2) = -cO * si;

        R(2, 0) =  sw * si;
        R(2, 1) =  cw * si;
        R(2, 2) =  ci;

        return R;
    }

    double OrbitPropagator::computeGMST() const {
        double gmst = gmst0_ + omega_earth_ * elapsed_time_;
        gmst = std::fmod(gmst, 2.0 * M_PI);
        if (gmst < 0.0) gmst += 2.0 * M_PI;
        return gmst;
    }

    Eigen::Vector3d OrbitPropagator::eciToECEF(const Eigen::Vector3d &r_eci) const {
        double gmst = computeGMST();
        double c = std::cos(gmst), s = std::sin(gmst);

        Eigen::Matrix3d Rz;
        Rz <<  c, s, 0,
              -s, c, 0,
               0, 0, 1;

        return Rz * r_eci;
    }

    LLA OrbitPropagator::ecefToLLA(const Eigen::Vector3d &r_ecef) const {
        double x = r_ecef.x(), y = r_ecef.y(), z = r_ecef.z();

        double a  = R_earth_;
        double e2 = 2.0 * f_ - f_ * f_;

        double lon = std::atan2(y, x);

        double p   = std::sqrt(x * x + y * y);
        double lat = std::atan2(z, p * (1.0 - e2));

        for (int i = 0; i < 10; ++i) {
            double sinLat = std::sin(lat);
            double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
            lat = std::atan2(z + e2 * N * sinLat, p);
        }

        double sinLat = std::sin(lat);
        double N   = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
        double alt = p / std::cos(lat) - N;

        LLA lla;
        lla.lat_deg = lat * 180.0 / M_PI;
        lla.lon_deg = lon * 180.0 / M_PI;
        lla.alt_km  = alt / 1000.0;
        return lla;
    }
}