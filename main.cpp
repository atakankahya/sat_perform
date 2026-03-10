#include <iostream>
#include <Eigen/Dense>
#include "spacecraft/SpacecraftState.h"

int main() {
    using namespace spacecraft;
    Eigen::Quaterniond q_init =
        Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(30, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(60, Eigen::Vector3d::UnitZ());
        Eigen::Vector3d omega_init(0.1,-0.1,0.1);

    SpacecraftState state(q_init, omega_init);
    std::cout << "States: " << state.DCM() << std::endl;
    std::cout << "Quat: " << state.q.norm() << std::endl;

    return 0;
}