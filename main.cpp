#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include "spacecraft/SpacecraftState.h"
#include "sensors/IGRF.h"
#include "orbit/OrbitPropagator.h"
#include "config/OrbitConfig.h"
#include "config/SimConfig.h"
#include "utils/PropagateDate.h"
#include "utils/GMST.h"
#include "sensors/Magnetosensor.h"
#include "spacecraft/Kinematics.h"
#include "utils/Conversions.h"
#include "config/SatelliteConfig.h"
#include "controller/BdotController.h"
#include <fstream>
#include <filesystem>
#include <stdexcept>

int main() {
    config::SimConfig sim;
    config::SatelliteConfig sat = config::SatelliteConfig::create(config::SatelliteType::UPMSAT2, false);
    double k = -3499767.38814494;
    double meas_interval = 0.2;

    controller::BdotController bdot(k, sat.max_dipole,meas_interval);
    config::OrbitConfig orb = config::OrbitConfig::create(600,78.7);
    const double gmst0 = utils::computeGMST(sim.start_date, 0.0);
    orbit::OrbitPropagator prop(orb, gmst0);

    std::filesystem::path sensors_dir = "sensors";
    const std::filesystem::path coeff_file = sensors_dir / "IGRF14coeffs.dat";
    if (!std::filesystem::exists(coeff_file)) {
        throw std::runtime_error(
            "IGRF coefficients not found at " + coeff_file.string() +
            ". Run the executable from the repository root (where sensors/ exists)."
        );
    }

    sensors::IGRF igrf(sensors_dir.string());
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


    std::ofstream csv("sim_log.csv");
    csv << "t,step_in_cycle,wx,wy,wz,w_norm,mx,my,mz,pwm_x,pwm_y,pwm_z,taux,tauy,tauz,"
           "Bx,By,Bz,B_norm,m_dot_B,tau_dot_omega,q0,q1,q2,q3\n";

    while (t < t_end) {
        orbit::LLA lla = prop.getLLA();
        double dec_year  = utils::advanceDecimalYear(decimal_year0,t);
        // state.DCM() is BODY->ECI with current quaternion convention.
        Eigen::Matrix3d R_eci2body = state.DCM().transpose();
        Eigen::Vector3d B_body = mag.measure(lla, prop.computeGMST(), dec_year, R_eci2body);


        int step_in_cycle = step % steps_per_cycle;

        // phase 1 measruement

        if (step_in_cycle < meas_phase_steps && step_in_cycle % meas_every == 0) {
            bdot.addMeasurement(B_body*1e-9);
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
        Eigen::Vector3d tau = controller::BdotController::torque(m_applied, B_body *1e-9);
        prop.propagate(dt);

        spacecraft::integrateOmega(state,sat.inertia,tau,dt);
        spacecraft::integrateEuler(state,dt);
        t += dt;

        if (1) {
            double m_dot_B = m_applied.dot(B_body * 1e-9);
            double tau_dot_omega = tau.dot(state.omega);
            csv << t << ","
                << step_in_cycle << ","
                << state.omega(0) << "," << state.omega(1) << "," << state.omega(2) << ","
                << state.omega.norm() << ","
                << m_applied(0) << "," << m_applied(1) << "," << m_applied(2) << ","
                << ctrl.pulse_ms(0) << "," << ctrl.pulse_ms(1) << "," << ctrl.pulse_ms(2) << ","
                << tau(0) << "," << tau(1) << "," << tau(2) << ","
                << B_body(0) << "," << B_body(1) << "," << B_body(2) << ","
                << (B_body * 1e-9).norm() << ","
                << m_dot_B << ","
                << tau_dot_omega << ","
                << state.q.w() << "," << state.q.x() << "," << state.q.y() << "," << state.q.z()
                << "\n";
        }

        step++;
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Initial omega: " << sim.omega_init.transpose()
              << "  norm = " << sim.omega_init.norm() << " rad/s\n";
    std::cout << "Final omega:   " << state.omega.transpose()
              << "  norm = " << state.omega.norm() << " rad/s\n";
    std::cout << "Reduction:     "
              << (1.0 - state.omega.norm() / sim.omega_init.norm()) * 100.0 << " %\n";

    csv.close();
    return 0;
}
