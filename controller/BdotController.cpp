//
// Created by ataka on 11.03.2026.
//

#include "BdotController.h"
#include "Eigen/Dense"
#include <cmath>


namespace controller {
    BdotController::BdotController(double gain, double max_dipole, double meas_interval_s)
        : k_(gain),
          max_dipole_(max_dipole),
          meas_interval_s_(meas_interval_s)
    {}

    void BdotController::addMeasurement(const Eigen::Vector3d &B_body) {
        measurements_.push_back(B_body);
    }

    ControlOutput BdotController::computeControl(const Eigen::Vector3d &omega_desired) {
        ControlOutput ctrl;
        ctrl.sign     = Eigen::Vector3d::Zero();
        ctrl.pulse_ms = Eigen::Vector3d::Zero();

        if (measurements_.size() < 2) return ctrl;
        // get 5 samples just like in munoz perez

        int n = static_cast<int> (measurements_.size());
        double total_time = (n-1) * meas_interval_s_;
        Eigen::Vector3d Bdot = (measurements_.back() - measurements_.front()) / total_time;

        //use last measruement as reference
        Eigen::Vector3d B_body = measurements_.back();

        Eigen::Vector3d signal = k_ * (omega_desired.cross(B_body)) + Bdot;

        for (int i = 0; i < 3; i++) {
            double mag = std::abs(signal(i));
            ctrl.pulse_ms(i) = pulseDuration(mag);

            if (signal(i) > 0) ctrl.sign(i) = 1;
            else if (signal(i) < 0) ctrl.sign(i) = -1;
            else  ctrl.sign(i) = 0;
        }

        return ctrl;

    }


    Eigen::Vector3d BdotController::getDipole(const ControlOutput &ctrl, double time_in_actuation_s) const {
        Eigen::Vector3d m = Eigen::Vector3d::Zero();
        for (int i = 0; i < 3; i++) {
            double pulse_s = ctrl.pulse_ms(i) / 1000.0;
            if (time_in_actuation_s < pulse_s) {
                m(i) = ctrl.sign(i) * max_dipole_;
            }
        }
        return m;
    }

    Eigen::Vector3d BdotController::torque(const Eigen::Vector3d& m_applied,
                                           const Eigen::Vector3d& B_body)
    {
        return m_applied.cross(B_body);
    }

    void BdotController::resetCycle() {
        measurements_.clear();
    }

    double BdotController::pulseDuration(double magnitude) {
        if (magnitude < 0.05) return 0.0;
        if (magnitude < 0.15) return 100.0;
        if (magnitude < 0.25) return 200.0;
        if (magnitude < 0.35) return 300.0;
        if (magnitude < 0.45) return 400.0;
        return 500.0;
    }
}
