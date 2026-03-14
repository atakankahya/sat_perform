//
// Created by ataka on 2.03.2026.
//

#ifndef SAT_PERFORM_CONVERSIONS_H
#define SAT_PERFORM_CONVERSIONS_H

#include "Eigen/Dense"
#include <cmath>

namespace utils {

    inline Eigen::Quaterniond eulerToQuaternion(double roll, double pitch, double yaw) {
        return Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX())
             * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
             * Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ());
    }

    inline Eigen::Quaterniond scalarFirst2VectorFirst(const Eigen::Vector4d& q) {
        return Eigen::Quaterniond(q(1), q(2), q(3), q(0));
    }

    inline Eigen::Quaterniond VectorFirst2ScalarFirst(Eigen::Vector4d& q) {
        return Eigen::Quaterniond(q(0), q(1), q(2), q(3));
    }

    /// Rotation matrix that maps a vector expressed in NED into ECEF.
    ///
    /// The columns are the NED unit-vectors written in ECEF:
    ///   col 0 = n_hat  (North)
    ///   col 1 = e_hat  (East)
    ///   col 2 = d_hat  (Down)
    ///
    /// NOTE: The original version returned the TRANSPOSE (ECEF→NED).
    ///       This version returns R_ned2ecef so that
    ///           v_ecef = R_ned2ecef * v_ned
    inline Eigen::Matrix3d nedToECEFDCM(double lat_deg, double lon_deg) {
        constexpr double D2R = M_PI / 180.0;
        double lat = lat_deg * D2R;
        double lon = lon_deg * D2R;
        double slat = std::sin(lat);
        double clat = std::cos(lat);
        double slon = std::sin(lon);
        double clon = std::cos(lon);

        // R_ecef2ned (rows = N, E, D expressed in ECEF)
        Eigen::Matrix3d R_ecef2ned;
        R_ecef2ned << -slat * clon,  -slat * slon,   clat,
                      -slon,          clon,           0.0,
                      -clat * clon,  -clat * slon,   -slat;

        // We need NED→ECEF, which is the transpose
        return R_ecef2ned.transpose();
    }

    /// Rotation matrix that maps a vector expressed in ECEF into ECI.
    ///   v_eci = R_ecef2eci * v_ecef
    /// This is Rz(+GMST).
    inline Eigen::Matrix3d ecefToEciDCM(double gmst_rad) {
        double c = std::cos(gmst_rad);
        double s = std::sin(gmst_rad);
        Eigen::Matrix3d R;
        R <<  c, -s,  0,
              s,  c,  0,
              0,  0,  1;
        return R;
    }

} // namespace utils

#endif // SAT_PERFORM_CONVERSIONS_H