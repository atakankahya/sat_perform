//
// Created by ataka on 10.03.2026.
// Faithful C++ port of MATLAB igrf.m
//

#include "IGRF.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sensors {

    IGRF::IGRF(const std::string& datadir)
    : baseEpoch_(2025.0)
{
    for (auto& row : g_)   for (auto& v : row) v = 0.0;
    for (auto& row : h_)   for (auto& v : row) v = 0.0;
    for (auto& row : gSV_) for (auto& v : row) v = 0.0;
    for (auto& row : hSV_) for (auto& v : row) v = 0.0;

    const std::filesystem::path path = std::filesystem::path(datadir) / "IGRF14coeffs.dat";
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("IGRF: cannot open " + path.string());

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string gh;
        int n, m;
        if (!(iss >> gh >> n >> m)) continue;
        if (gh != "g" && gh != "h") continue;
        if (n < 1 || n > NMAX || m < 0 || m > n) continue;

        std::vector<double> vals;
        double v;
        while (iss >> v) vals.push_back(v);
        if (vals.size() < 2) continue;

        double main_val = vals[vals.size() - 2];
        double sv_val   = vals[vals.size() - 1];

        if (gh == "g") {
            g_[n][m]   = main_val;
            gSV_[n][m] = sv_val;
        } else {
            h_[n][m]   = main_val;
            hSV_[n][m] = sv_val;
        }
    }
}

// ── computeNED — faithful port of MATLAB igrf.m ──────────────────────────
Eigen::Vector3d IGRF::computeNED(double lat_deg, double lon_deg,
                                  double alt_km, double decimal_year) const
{
    constexpr double D2R = M_PI / 180.0;
    constexpr double R2D = 180.0 / M_PI;
    const int nmax = NMAX;  // 13

    const double dt = decimal_year - baseEpoch_;
    double gc[N][N]{}, hc[N][N]{};
    for (int n = 1; n <= nmax; ++n)
        for (int m = 0; m <= n; ++m) {
            gc[n][m] = g_[n][m] + gSV_[n][m] * dt;
            hc[n][m] = h_[n][m] + hSV_[n][m] * dt;
        }

    const double b_axis = WGS84_a * std::sqrt(1.0 - WGS84_e2);  // semi-minor axis

    const double gdlat_rad = lat_deg * D2R;
    const double sin_alpha_2 = std::sin(gdlat_rad) * std::sin(gdlat_rad);
    const double cos_alpha_2 = std::cos(gdlat_rad) * std::cos(gdlat_rad);

    const double tmp = alt_km * std::sqrt(WGS84_a * WGS84_a * cos_alpha_2
                                        + b_axis * b_axis * sin_alpha_2);

    const double beta  = std::atan2((tmp + b_axis * b_axis) * std::tan(gdlat_rad),
                                     tmp + WGS84_a * WGS84_a);
    const double theta = M_PI / 2.0 - beta;               // geocentric colatitude

    const double ba2 = (b_axis / WGS84_a) * (b_axis / WGS84_a);
    const double ba4 = ba2 * ba2;
    const double r_km = std::sqrt(
        alt_km * alt_km + 2.0 * tmp
        + WGS84_a * WGS84_a * (1.0 - (1.0 - ba4) * sin_alpha_2)
                              / (1.0 - (1.0 - ba2) * sin_alpha_2));

    const double theta_deg = theta * R2D;
    const double phi_deg   = lon_deg;

    const double theta_rad = theta;  // already in radians
    const double costh = std::cos(theta_rad);
    const double sinth = std::sin(theta_rad);

    double P [N][N]{};
    double dP[N][N]{};
    double S [N][N]{};

    P [0][0] = 1.0;
    dP[0][0] = 0.0;
    S [0][0] = 1.0;

    for (int n = 1; n <= nmax; ++n) {
        for (int m = 0; m <= n; ++m) {
            if (n == m) {
                // diagonal
                P [n][m] = sinth * P[n - 1][m - 1];
                dP[n][m] = sinth * dP[n - 1][m - 1] + costh * P[n - 1][n - 1];
            } else {
                // off-diagonal
                double Knm = 0.0;
                if (n > 1) {
                    Knm = (static_cast<double>((n - 1) * (n - 1) - m * m))
                        / (static_cast<double>((2 * n - 1) * (2 * n - 3)));
                }

                P [n][m] = costh * P[n - 1][m];
                dP[n][m] = costh * dP[n - 1][m] - sinth * P[n - 1][m];

                if (n > 1) {
                    P [n][m] -= Knm * P [n - 2][m];
                    dP[n][m] -= Knm * dP[n - 2][m];
                }
            }

            // Schmidt normalisation factor S
            if (m == 0) {
                S[n][0] = S[n - 1][0] * (2.0 * n - 1.0) / n;
            } else {
                double extra = (m == 1) ? 2.0 : 1.0;
                S[n][m] = S[n][m - 1]
                         * std::sqrt((n - m + 1.0) * extra / (n + m));
            }
        }
    }

    // Apply Schmidt factors
    for (int n = 0; n <= nmax; ++n)
        for (int m = 0; m <= n; ++m) {
            P [n][m] *= S[n][m];
            dP[n][m] *= S[n][m];
        }

    double sinth_safe = sinth;
    if (std::abs(sinth_safe) < 1.0e-12)
        sinth_safe = (sinth_safe >= 0.0 ? 1.0e-12 : -1.0e-12);

    double Br = 0.0, Btheta = 0.0, Bphi = 0.0;

    for (int n = 1; n <= nmax; ++n) {
        const double ratio_n = std::pow(RE / r_km, n + 2);

        for (int m = 0; m <= n; ++m) {
            const double gnm = gc[n][m];
            const double hnm = hc[n][m];

            const double cm = std::cos(m * phi_deg * D2R);
            const double sm = std::sin(m * phi_deg * D2R);

            const double C   = gnm * cm + hnm * sm;
            const double Pnm  = P [n][m];
            const double dPnm = dP[n][m];

            Br     += ratio_n * (n + 1) * C * Pnm;
            Btheta -= ratio_n * C * dPnm;

            if (m > 0) {
                Bphi += ratio_n * m * (gnm * sm - hnm * cm) * Pnm / sinth_safe;
            }
        }
    }

    const double E2 = 1.0 - ba2;
    const double E4 = E2 * E2;
    const double E6 = E4 * E2;
    const double E8 = E4 * E4;

    const double A21 =  (512*E2 + 128*E4 + 60*E6 + 35*E8) / 1024.0;
    const double A22 =  (E6 + E8) / 32.0;
    const double A23 = -3.0 * (4*E6 + 3*E8) / 256.0;
    const double A41 = -(64*E4 + 48*E6 + 35*E8) / 1024.0;
    const double A42 =  (4*E4 + 2*E6 + E8) / 16.0;
    const double A43 =  15.0 * E8 / 256.0;
    const double A44 = -E8 / 16.0;
    const double A61 =  3.0 * (4*E6 + 5*E8) / 1024.0;
    const double A62 = -3.0 * (E6 + E8) / 32.0;
    const double A63 =  35.0 * (4*E6 + 3*E8) / 768.0;
    const double A81 = -5.0 * E8 / 2048.0;
    const double A82 =  64.0 * E8 / 2048.0;
    const double A83 = -252.0 * E8 / 2048.0;
    const double A84 =  320.0 * E8 / 2048.0;

    const double GCLAT = 90.0 - theta_deg;
    const double SCL   = std::sin(GCLAT * D2R);

    const double RI = WGS84_a / r_km;
    const double A2 = RI * (A21 + RI * (A22 + RI * A23));
    const double A4 = RI * (A41 + RI * (A42 + RI * (A43 + RI * A44)));
    const double A6 = RI * (A61 + RI * (A62 + RI * A63));
    const double A8 = RI * (A81 + RI * (A82 + RI * (A83 + RI * A84)));

    const double CCL  = std::sqrt(1.0 - SCL * SCL);
    const double S2CL = 2.0 * SCL * CCL;
    const double C2CL = 2.0 * CCL * CCL - 1.0;
    const double S4CL = 2.0 * S2CL * C2CL;
    const double C4CL = 2.0 * C2CL * C2CL - 1.0;
    const double S8CL = 2.0 * S4CL * C4CL;
    const double S6CL = S2CL * C4CL + C2CL * S4CL;

    const double DLTCL     = S2CL * A2 + S4CL * A4 + S6CL * A6 + S8CL * A8;
    const double gdlat_calc = GCLAT * D2R + DLTCL;

    const double psi = std::sin(gdlat_calc) * std::sin(theta_rad)
                     - std::cos(gdlat_calc) * std::cos(theta_rad);

    const double Bn_temp = -std::cos(psi) * Btheta - std::sin(psi) * Br;
    const double Bu_temp = -std::sin(psi) * Btheta + std::cos(psi) * Br;
    const double Be_temp =  Bphi;

    // NED: North, East, Down (Down = -Up)
    const double BN =  Bn_temp;
    const double BE =  Be_temp;
    const double BD = -Bu_temp;

    return {BN, BE, BD};
}

} // namespace sensors
