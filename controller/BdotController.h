//
// Created by ataka on 11.03.2026.
//

#ifndef SAT_PERFORM_BDOTCONTROLLER_H
#define SAT_PERFORM_BDOTCONTROLLER_H
#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include <vector>

namespace controller {

    struct ControlOutput {
        Eigen::Vector3d sign;
        Eigen::Vector3d pulse_ms;
    };


    class BdotController {
    public:
        BdotController(double gain, double max_dipole, double meas_interval_s);


        void addMeasurement(const Eigen::Vector3d& B_body);


        // 5 measruement and compute the sign and pwm per axis.
        ControlOutput computeControl(const Eigen::Vector3d& omega_desired);
        Eigen::Vector3d getDipole(const ControlOutput& ctrl, double time_in_actuation_s) const;

        static Eigen::Vector3d torque(const Eigen::Vector3d& m_applied,
                                      const Eigen::Vector3d& B_body);

        void resetCycle();


    private:
        double k_;
        double max_dipole_;
        double meas_interval_s_;
        std::vector<Eigen::Vector3d> measurements_;


        static double pulseDuration(double magnitude);

    };

}


#endif //SAT_PERFORM_BDOTCONTROLLER_H