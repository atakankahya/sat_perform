//
// Created by ataka on 14.03.2026.
//

#ifndef SAT_PERFORM_SATENV_H
#define SAT_PERFORM_SATENV_H


#include <array>
#include "Eigen/Dense"
#include <random>
#include <vector>
#include <deque>


#include "spacecraft/SpacecraftState.h"
#include "sensors/IGRF.h"
#include "sensors/Magnetosensor.h"
#include "config/OrbitConfig.h"
#include "config/SimConfig.h"
#include "config/SatelliteConfig.h"
#include "utils/Conversions.h"
#include "utils/PropagateDate.h"
#include "utils/GMST.h"
#include "spacecraft/Kinematics.h"




class SatEnv {
public:
    static constexpr int OBS_DIM = 9;
    static constexpr int N_ACTIONS = 125;
    static constexpr int LEVELS_PER_AXIS = 5;

    explicit SatEnv(const std::string& igrf_datadir);

    std::array<double, OBS_DIM> reset();
    std::tuple<std::array<double, OBS_DIM>, double, bool> step(int action);

private:
    config::SimConfig sim_;
    config::OrbitConfig orb_;
    config::SatelliteConfig sat_;

    sensors::IGRF igrf_;
    sensors::MagnetoSensor mag_;


    std::unique_ptr<orbit::OrbitPropagator> prop_;
    spacecraft::SpacecraftState state_;

    double t_;
    double decimal_year0_;
    double gmst0_;
    int episode_steps_;
    int max_episode_steps_;


    std::mt19937 rng_;

    static constexpr double PWM_LEVELS[5] = {-500.0, -200.0, 0.0, 200.0, 500.0};
    static constexpr double MEAS_INTERVAL_S = 0.2;
    static constexpr double MEAS_STEPS = 5;


    static constexpr int REWARD_WINDOW = 50.0;
    std::deque<double> reward_history_;


    //minmaxtable from munoz paper
    static constexpr double BDOT_MIN = -8E6;
    static constexpr double BDOT_MAX = 8E6;
    static constexpr double BMEAN_MIN = -5E5;
    static constexpr double BMEAN_MAX = 5E5;
    static constexpr double WERR_MIN = 0.0;
    static constexpr double WERR_MAX  = 0.25;

    void decodeAction(int action, Eigen::Vector3d& pwm_ms) const;
    std::array<double, OBS_DIM> buildObs(const Eigen::Vector3d& Bdot_nT,const Eigen::Vector3d& Bmean_nT) const;

    double computeBaseReward() const;
    double computeFilteredReward(double r_base);
    Eigen::Vector3d measureB_nT();

    static double normalize(double x, double xmin, double xmax);










public:




};


#endif //SAT_PERFORM_SATENV_H