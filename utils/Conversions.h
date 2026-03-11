//
// Created by ataka on 2.03.2026.
//

#ifndef SAT_PERFORM_CONVERSIONS_H
#define SAT_PERFORM_CONVERSIONS_H

#include "Eigen/Dense"

namespace utils {
    inline Eigen::Quaterniond eulerToQuaternion(double roll, double pitch, double yaw) {
        return Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX())
             * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
             * Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ());
    }

    inline Eigen::Quaterniond scalarFirst2VectorFirst(const Eigen::Vector4d& q) {
        return Eigen::Quaterniond(q(1),q(2),q(3),q(0));
    }

    inline Eigen::Quaterniond VectorFirst2ScalarFirst(Eigen::Vector4d& q) {
        return Eigen::Quaterniond(q(0),q(1),q(2),q(3));
    }

    //NED to ecef dcm
    inline Eigen::Matrix3d nedToECEFDCM(double lat_deg, double lon_deg) {
        constexpr double D2R = M_PI / 180.0;
        double lat = lat_deg * D2R;
        double lon = lon_deg * D2R;
        double slat = std::sin(lat);
        double clat = std::cos(lat);
        double slon = std::sin(lon);
        double clon = std::cos(lon);

        Eigen::Matrix3d R;
        R << -slat * clon,   -slon,   -clat * clon,
        -slat * slon,    clon,   -clat * slon,
         clat,            0.0,   -slat;
        return R;
    }

    //ecef to eci rotation matrix
    inline Eigen::Matrix3d ecefToEciDCM(double gmst_rad) {
        double c = std::cos(gmst_rad);
        double s = std::sin(gmst_rad);
        Eigen::Matrix3d R;
        R <<  c, -s,  0,
        s,  c,  0,
        0,  0,  1;

        return R;

    }

}


#endif //SAT_PERFORM_CONVERSIONS_H