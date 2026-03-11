#include <iostream>
#include <Eigen/Dense>
#include "spacecraft/SpacecraftState.h"
#include "sensors/IGRF.h"
#include "orbit/OrbitPropagator.h"
#include "config/OrbitConfig.h"
#include "utils/PropagateDate.h"


int main() {

    sensors::IGRF igrf("sensors");
    double decimal_year = utils::getDecimalYear();


    Eigen::Vector3d B_ground = igrf.computeNED(-23.5, -46.6, 770.0, decimal_year);
    std::cout << "B ground (nT):  " << B_ground.transpose() << std::endl;
    std::cout << "|B| = " << B_ground.norm() << " nT\n" << std::endl;




    return 0;
}
