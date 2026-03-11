#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include "spacecraft/SpacecraftState.h"
#include "sensors/IGRF.h"
#include "orbit/OrbitPropagator.h"
#include "config/OrbitConfig.h"
#include "config/SimConfig.h"
#include "utils/PropagateDate.h"
#include "sensors/Magnetosensor.h"
#include "spacecraft/Kinematics.h"
#include "utils/Conversions.h"


int main() {
    config::SimConfig sim;
    config::OrbitConfig orb = config::OrbitConfig::create(550,97.4);
    orbit::OrbitPropagator prop(orb);

    sensors::IGRF igrf("sensors");
    sensors::MagnetoSensor mag(igrf);

    Eigen::Quaterniond q0 = utils::eulerToQuaternion(
        sim.roll_deg * M_PI/180,
        sim.pitch_deg * M_PI/180,
        sim.yaw_deg * M_PI/180);

    // Change the initial conditions from here.

    spacecraft::SpacecraftState state(q0, sim.omega_init);

    double decimal_year0 = utils::dateToDecimalYear({
    static_cast<int>(sim.start_date[0]),
    static_cast<int>(sim.start_date[1]),
    static_cast<int>(sim.start_date[2]),
    static_cast<int>(sim.start_date[3]),
    static_cast<int>(sim.start_date[4]),
    sim.start_date[5]
    });

    double dt = sim.dt;
    double t_end = sim.t_episode;
    double t = 0.0;
    int step = 0;
    int print_every = static_cast<int>(sim.control_period / dt);

    std::cout << std::fixed << std::setprecision(2) << std::setw(10) << t << std::endl;
    std::cout << "Time (s) " << std::endl;


    while (t < t_end) {
        orbit::LLA lla = prop.getLLA();
        double dec_year  = utils::advanceDecimalYear(decimal_year0,t);
        Eigen::Matrix3d Rb = state.DCM();

        Eigen::Vector3d B_body = mag.measure(lla, prop.computeGMST(),dec_year, Rb);


        if (step % print_every == 0) {
            std::cout << std::setw(8) << t << std::endl;
            std::cout << std::setw(12) << B_body(0) << std::endl;
            std::cout << std::setw(12) << B_body(1) << std::endl;
            std::cout << std::setw(12) << B_body(2) << std::endl;
        }
        prop.propagate(dt);

        spacecraft::integrateEuler(state,dt);
        t += dt;
        step++;

    }

    return 0;
}
