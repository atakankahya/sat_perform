//
// Created by ataka on 10.03.2026.
//

#include "IGRF.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace sensors {

    IGRF::IGRF(const std::string &datadir)
        : baseEpoch_(2025.0)
    {
        for (auto& row : g_)   for (auto& v : row) v = 0.0;
        for (auto& row : h_)   for (auto& v : row) v = 0.0;
        for (auto& row : gSV_) for (auto& v : row) v = 0.0;
        for (auto& row : hSV_) for (auto& v : row) v = 0.0;

        const std::string path = datadir + "/IGRF14coeffs.dat";
        std::ifstream file(path);
        if (!file.is_open())
            throw std::runtime_error("IGRF: cannot open " + path);

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            std::string gh;
            int n, m;
            if (!(iss >> gh >> n >> m)) continue;
            if (gh != "g" && gh != "h") continue;
            if (n < 1 || n > NMAX || m < 0 || m > n) continue;

            // read all remaining doubles — second-to-last is 2025.0, last is SV
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


    Eigen::Vector3d IGRF::computeNED(double lat_deg, double lon_deg,
                                      double alt_km, double decimal_year) const
    {
        // interpolate coefficients to requested epoch
        const double dt = decimal_year - baseEpoch_;
        double gc[N][N]{}, hc[N][N]{};
        for (int n = 1; n <= NMAX; ++n)
            for (int m = 0; m <= n; ++m) {
                gc[n][m] = g_[n][m] + gSV_[n][m] * dt;
                hc[n][m] = h_[n][m] + hSV_[n][m] * dt;
            }

        // geodetic to geocentric (WGS-84)
        constexpr double D2R = M_PI / 180.0;
        const double lat = lat_deg * D2R;
        const double lon = lon_deg * D2R;
        const double slat = std::sin(lat);
        const double clat = std::cos(lat);

        const double Rn  = WGS84_a / std::sqrt(1.0 - WGS84_e2 * slat * slat);
        const double rho = (Rn + alt_km) * clat;
        const double z   = (Rn * (1.0 - WGS84_e2) + alt_km) * slat;
        const double r   = std::sqrt(rho * rho + z * z);
        const double ct  = z / r;
        const double st  = rho / r;

        const double gc_colat = std::atan2(rho, z);
        const double gd_colat = M_PI / 2.0 - lat;
        const double delta    = gd_colat - gc_colat;
        const double cd = std::cos(delta);
        const double sd = std::sin(delta);

        // Schmidt semi-normalised Legendre polynomials and dP/dtheta
        double P [N][N]{};
        double dP[N][N]{};

        P [0][0] = 1.0;   dP[0][0] =  0.0;
        P [1][0] = ct;     dP[1][0] = -st;
        P [1][1] = st;     dP[1][1] =  ct;

        for (int n = 2; n <= NMAX; ++n) {
            // diagonal m == n
            double dn = std::sqrt((2.0 * n - 1.0) / (2.0 * n));
            P [n][n] = dn * st * P[n-1][n-1];
            dP[n][n] = dn * (ct * P[n-1][n-1] + st * dP[n-1][n-1]);

            // sub-diagonal m == n-1
            double sn = std::sqrt(2.0 * n - 1.0);
            P [n][n-1] = sn * ct * P[n-1][n-1];
            dP[n][n-1] = sn * (-st * P[n-1][n-1] + ct * dP[n-1][n-1]);

            // general 0 <= m <= n-2
            for (int m = 0; m <= n - 2; ++m) {
                double nd  = static_cast<double>(n);
                double md  = static_cast<double>(m);
                double den = std::sqrt(nd*nd - md*md);
                double a   = (2.0*nd - 1.0) / den;
                double b   = std::sqrt((nd-1.0)*(nd-1.0) - md*md) / den;

                P [n][m] = a * ct * P[n-1][m] - b * P[n-2][m];
                dP[n][m] = a * (-st * P[n-1][m] + ct * dP[n-1][m]) - b * dP[n-2][m];
            }
        }

        // spherical harmonic summation
        double Br = 0.0, Bt = 0.0, Bp = 0.0;
        double rn = (RE / r) * (RE / r);

        for (int n = 1; n <= NMAX; ++n) {
            rn *= (RE / r);
            for (int m = 0; m <= n; ++m) {
                double cm = std::cos(m * lon);
                double sm = std::sin(m * lon);
                double T  = gc[n][m] * cm + hc[n][m] * sm;
                double Q  = m * (-gc[n][m] * sm + hc[n][m] * cm);

                Br += static_cast<double>(n + 1) * rn * T * P[n][m];
                Bt -= rn * T * dP[n][m];
                if (std::abs(st) > 1.0e-10)
                    Bp -= rn * Q * P[n][m] / st;
            }
        }

        // geocentric spherical to geodetic NED
        double BN = -Bt * cd - Br * sd;
        double BE =  Bp;
        double BD =  Bt * sd - Br * cd;

        return {BN, BE, BD};
    }

}
