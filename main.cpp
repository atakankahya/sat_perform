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
#include "config/SatelliteConfig.h"
#include "controller/BdotController.h"


int main() {
    config::SimConfig sim;
    config::SatelliteConfig sat = config::SatelliteConfig::create(config::SatelliteType::Cubesat3U, false);
    double k = 3499767.38814494;
    double meas_interval = 0.2;

    controller::BdotController bdot(k, sat.max_dipole,meas_interval);
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
    double meas_phase = 1.0;
    double control_period = sim.control_period;

    int steps_per_cycle  = static_cast<int>(control_period  / dt);   // 20
    int meas_phase_steps = static_cast<int>(meas_phase    / dt);   // 10
    int meas_every       = static_cast<int>(meas_interval / dt);   //  2

    controller::ControlOutput ctrl;
    ctrl.sign     = Eigen::Vector3d::Zero();
    ctrl.pulse_ms = Eigen::Vector3d::Zero();


    while (t < t_end) {
        orbit::LLA lla = prop.getLLA();
        double dec_year  = utils::advanceDecimalYear(decimal_year0,t);
        Eigen::Matrix3d Rb = state.DCM();
        Eigen::Vector3d B_body = mag.measure(lla, prop.computeGMST(),dec_year, Rb);


        int step_in_cycle = step % steps_per_cycle;

        // phase 1 measruement

        if (step_in_cycle < meas_phase_steps && step_in_cycle % meas_every == 0) {            bdot.addMeasurement(B_body);
        }

        //control phase 2
        if (step_in_cycle == meas_phase_steps) {
            ctrl = bdot.computeControl(sim.omega_desired);
            bdot.resetCycle();
        }

        // actuation phase
        Eigen::Vector3d m_applied = Eigen::Vector3d::Zero();
        if (step_in_cycle >= meas_phase_steps) {
            double time_in_actuation = (step_in_cycle - meas_phase_steps) * dt;
            m_applied = bdot.getDipole(ctrl,time_in_actuation);
        }

        //dynamics
        Eigen::Vector3d tau = controller::BdotController::torque(m_applied, B_body);
        prop.propagate(dt);

        spacecraft::integrateOmega(state,sat.inertia,tau,dt);
        spacecraft::integrateEuler(state,dt);
        t += dt;
        step++;

    }

    return 0;
}
