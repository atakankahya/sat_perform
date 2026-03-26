//
// Created by ataka on 14.03.2026.
//

#include "SatEnv.h"
#include <algorithm>

SatEnv::SatEnv(const std::string &igrf_datadir)
    :  sat_(config::SatelliteConfig::create(config::SatelliteType::UPMSAT2, false))
    ,  orb_(config::OrbitConfig::create(600, 78.1))
    ,  igrf_(igrf_datadir)
    ,  mag_(igrf_)
    ,  t_(0.0)
    ,  decimal_year0_(0.0)
    ,  gmst0_(0.0)
    ,  episode_steps_(0.0)
    ,  max_episode_steps_(static_cast<int>(sim_.t_episode / sim_.control_period))
    ,  rng_(std::random_device{}())
{
}


double SatEnv::normalize(double x, double xmin, double xmax) {
    double val = (x- xmin) / (xmax - xmin);
    return std::clamp(val,0.0,1.0);
}


Eigen::Vector3d SatEnv::measureB_nT() {
    orbit::LLA lla = prop_->getLLA();
    double dec_year = utils::advanceDecimalYear(decimal_year0_,t_);
    Eigen::Matrix3d R_eci2body = state_.DCM().transpose();
    return mag_.measure(lla, prop_->computeGMST(), dec_year, R_eci2body);
}


void SatEnv::decodeAction(int action, Eigen::Vector3d& pwm_ms) const {
    // action 0..124 -> 3 indices into PWM_LEVELS {-500,-200,0,+200,+500}
    int a2 = action % LEVELS_PER_AXIS;
    int a1 = (action / LEVELS_PER_AXIS) % LEVELS_PER_AXIS;
    int a0 = action / (LEVELS_PER_AXIS * LEVELS_PER_AXIS);
    pwm_ms(0) = PWM_LEVELS[a0];
    pwm_ms(1) = PWM_LEVELS[a1];
    pwm_ms(2) = PWM_LEVELS[a2];
}

double SatEnv::computeBaseReward() const {
    // Eq (8)-(10): alpha = sum of absolute angular velocity errors
    Eigen::Vector3d w_err = sim_.omega_desired - state_.omega;
    double alpha = std::abs(w_err(0)) + std::abs(w_err(1)) + std::abs(w_err(2));
    return 1.0 - std::exp(2.0 * M_PI * alpha);
}

double SatEnv::computeFilteredReward(double r_base) {
    // Eq (11): moving average filter, window = 100s = 50 control steps
    reward_history_.push_back(r_base);
    if (static_cast<int>(reward_history_.size()) > REWARD_WINDOW)
        reward_history_.pop_front();

    double sum = 0.0;
    for (double r : reward_history_) sum += r;
    return sum / static_cast<double>(reward_history_.size());
}

std::array<double, SatEnv::OBS_DIM> SatEnv::buildObs(
    const Eigen::Vector3d& Bdot_nT,
    const Eigen::Vector3d& Bmean_nT) const
{
    // Table 3: 9 inputs, all min-max normalized to [0,1]
    Eigen::Vector3d w_err = sim_.omega_desired - state_.omega;

    return {
        normalize(Bdot_nT(0),  BDOT_MIN,  BDOT_MAX),
        normalize(Bdot_nT(1),  BDOT_MIN,  BDOT_MAX),
        normalize(Bdot_nT(2),  BDOT_MIN,  BDOT_MAX),
        normalize(Bmean_nT(0), BMEAN_MIN, BMEAN_MAX),
        normalize(Bmean_nT(1), BMEAN_MIN, BMEAN_MAX),
        normalize(Bmean_nT(2), BMEAN_MIN, BMEAN_MAX),
        normalize(w_err(0), -WERR_MAX, WERR_MAX),
        normalize(w_err(1), -WERR_MAX, WERR_MAX),
        normalize(w_err(2), -WERR_MAX, WERR_MAX),
            };
}



std::array<double, SatEnv::OBS_DIM> SatEnv::reset() {

    Eigen::Vector3d w0(0.1, -0.1, -0.1);
    Eigen::Quaterniond q0 = utils::eulerToQuaternion(
        0.0, 30.0 * M_PI / 180.0, 60.0 * M_PI / 180.0);

    state_ = spacecraft::SpacecraftState(q0, w0);

    // randomize orbital position
    std::uniform_real_distribution<double> nu_dist(0.0, 360.0);
    orb_.nu_deg = nu_dist(rng_);

    gmst0_ = utils::computeGMST(sim_.start_date, 0.0);
    prop_ = std::make_unique<orbit::OrbitPropagator>(orb_, gmst0_);

    decimal_year0_ = utils::dateToDecimalYear({
        static_cast<int>(sim_.start_date[0]),
        static_cast<int>(sim_.start_date[1]),
        static_cast<int>(sim_.start_date[2]),
        static_cast<int>(sim_.start_date[3]),
        static_cast<int>(sim_.start_date[4]),
        sim_.start_date[5]
    });

    t_ = 0.0;
    episode_steps_ = 0;
    reward_history_.clear();

    // return initial observation with zero Bdot and a single B measurement
    Eigen::Vector3d B0 = measureB_nT();
    return buildObs(Eigen::Vector3d::Zero(), B0);
}

// step

std::tuple<std::array<double, SatEnv::OBS_DIM>, double, bool> SatEnv::step(int action) {
    Eigen::Vector3d pwm_ms;
    decodeAction(action, pwm_ms);

    double dt = sim_.dt;
    int steps_per_cycle = static_cast<int>(sim_.control_period / dt); // 20

    // ── Phase 1: Measurement (first 1s = 10 substeps) ──
    // collect 5 magnetometer readings at 0.0, 0.2, 0.4, 0.6, 0.8 s
    int meas_phase_substeps = static_cast<int>(1.0 / dt); // 10
    int meas_every = static_cast<int>(MEAS_INTERVAL_S / dt); // 2

    std::vector<Eigen::Vector3d> B_samples;
    B_samples.reserve(MEAS_STEPS);

    for (int s = 0; s < meas_phase_substeps; s++) {
        if (s % meas_every == 0) {
            B_samples.push_back(measureB_nT());
        }
        // no actuation during measurement phase — just propagate
        prop_->propagate(dt);
        spacecraft::integrateOmega(state_, sat_.inertia, Eigen::Vector3d::Zero(), dt);
        spacecraft::integrateEuler(state_, dt);
        t_ += dt;
    }

    // compute Bmean and Bdot from the 5 samples
    Eigen::Vector3d Bmean = Eigen::Vector3d::Zero();
    for (auto& b : B_samples) Bmean += b;
    Bmean /= static_cast<double>(B_samples.size());

    // Bdot via linear regression (same as your BdotController)
    int N = static_cast<int>(B_samples.size());
    double Sx = 0, Sxx = 0;
    Eigen::Vector3d Sxy = Eigen::Vector3d::Zero();
    Eigen::Vector3d Sy  = Eigen::Vector3d::Zero();
    for (int i = 0; i < N; i++) {
        double x = i * MEAS_INTERVAL_S;
        Sx  += x;
        Sxx += x * x;
        Sxy += x * B_samples[i];
        Sy  += B_samples[i];
    }
    double denom = N * Sxx - Sx * Sx;
    Eigen::Vector3d Bdot = Eigen::Vector3d::Zero();
    if (std::abs(denom) > 1e-12) {
        Bdot = (N * Sxy - Sx * Sy) / denom;
    }

    // ── Phase 2: Actuation (remaining substeps) ──
    int actuation_substeps = steps_per_cycle - meas_phase_substeps;
    for (int s = 0; s < actuation_substeps; s++) {
        Eigen::Vector3d B_body_T = measureB_nT() * 1e-9;

        Eigen::Vector3d m_applied = Eigen::Vector3d::Zero();
        double time_in_actuation_s = s * dt;
        for (int i = 0; i < 3; i++) {
            double pulse_s = std::abs(pwm_ms(i)) / 1000.0;
            if (time_in_actuation_s < pulse_s) {
                double sign = (pwm_ms(i) > 0) ? 1.0 : ((pwm_ms(i) < 0) ? -1.0 : 0.0);
                m_applied(i) = sign * sat_.max_dipole;
            }
        }

        Eigen::Vector3d tau = m_applied.cross(B_body_T);
        prop_->propagate(dt);
        spacecraft::integrateOmega(state_, sat_.inertia, tau, dt);
        spacecraft::integrateEuler(state_, dt);
        t_ += dt;
    }

    episode_steps_++;
    bool done = (episode_steps_ >= max_episode_steps_);

    double r_base = computeBaseReward();
    double reward = computeFilteredReward(r_base);





    return {buildObs(Bdot, Bmean), reward, done};
}
