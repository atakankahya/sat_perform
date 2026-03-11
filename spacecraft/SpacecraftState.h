//
// Created by ataka on 27.02.2026.
//

#ifndef SAT_PERFORM_SPACECRAFTSTATE_H
#define SAT_PERFORM_SPACECRAFTSTATE_H

#include "Eigen/Dense"



namespace spacecraft {
    struct SpacecraftState {
        Eigen::Quaterniond q; // Quaternion and omega speed.
        Eigen::Vector3d omega;

        SpacecraftState()
            : q(Eigen::Quaterniond::Identity())
            , omega(Eigen::Vector3d::Zero())
        {}

        SpacecraftState(Eigen::Quaterniond q, Eigen::Vector3d omega)
            :q(std::move(q))
            , omega(std::move(omega))
        {}


        void normalizeQuaternion() {
            q.normalize();
        }

        // Quaternion to DCM conversion b_eci to b_body
        // i dont need to build something like in matlab calling eigen enough
        Eigen::Matrix3d DCM() const {
            return q.toRotationMatrix();
        }



    };
}


#endif //SAT_PERFORM_SPACECRAFTSTATE_H