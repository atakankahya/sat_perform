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

}


#endif //SAT_PERFORM_CONVERSIONS_H