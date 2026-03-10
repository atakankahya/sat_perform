//
// Created by ataka on 2.03.2026.
//

#include "Kinematics.h"
#include "SpacecraftState.h"


namespace spacecraft {
    Eigen::Matrix4d OmegaMatrix(const Eigen::Vector3d &omega) {

        double wx = omega(0);
        double wy = omega(1);
        double wz = omega(2);

        Eigen::Matrix4d Omega;

        Omega << 0, wz, -wy, wx,
        -wz, 0, wx, wy,
        wy, -wx, 0, wz,
        -wx, -wy, -wz, 0;

        return Omega;

    }


    Eigen::Vector4d computeQdot(const spacecraft::SpacecraftState &state) {
        Eigen::Vector4d q_vec = state.q.coeffs();
        Eigen::Matrix4d Omega = OmegaMatrix(state.omega);

        Eigen::Vector4d q_dot = 0.5 * Omega * q_vec;

        return q_dot;
    }


    void integrateEuler(spacecraft::SpacecraftState &state, double dt) {
        Eigen::Vector4d q_dot = computeQdot(state);
        Eigen::Vector4d q_vec = state.q.coeffs();
        Eigen::Vector4d q_new = q_vec + q_dot * dt;
        q_new.normalize();
        state.q = Eigen::Quaterniond(q_new(3),q_new(0),q_new(1),q_new(2));
    }
}