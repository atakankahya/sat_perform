//
// Created by ataka on 2.03.2026.
//

#ifndef SAT_PERFORM_KINEMATICS_H
#define SAT_PERFORM_KINEMATICS_H

#include "SpacecraftState.h"
#include "Eigen/Dense"




static Eigen::Matrix4d OmegaMatrix(const Eigen::Vector3d& omega);

static Eigen::Vector4d computeQdot(const spacecraft::SpacecraftState& state );

void integrateEuler(spacecraft::SpacecraftState& state, double dt);


#endif //SAT_PERFORM_KINEMATICS_H
