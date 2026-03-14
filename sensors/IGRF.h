//
// Created by ataka on 10.03.2026.
// Fixed constants and faithful port from MATLAB igrf.m
//

#ifndef SAT_PERFORM_IGRF_H
#define SAT_PERFORM_IGRF_H

#include "Eigen/Dense"
#include <string>

namespace sensors {

    class IGRF {
    public:
        explicit IGRF(const std::string& datadir);

        /// Compute magnetic field in NED frame (nanoTesla)
        /// @param lat_deg  geodetic latitude  [deg]
        /// @param lon_deg  geodetic longitude [deg]
        /// @param alt_km   altitude above WGS-84 ellipsoid [km]
        /// @param decimal_year  e.g. 2025.23
        [[nodiscard]] Eigen::Vector3d computeNED(double lat_deg,
                                                 double lon_deg,
                                                 double alt_km,
                                                 double decimal_year) const;

    private:
        static constexpr int NMAX = 13;
        static constexpr int N    = NMAX + 1;   // 14

        // IGRF reference radius
        static constexpr double RE       = 6371.2;            // km
        // WGS-84 ellipsoid
        static constexpr double WGS84_a  = 6378.137;          // km
        static constexpr double WGS84_e2 = 0.00669437999014;  // eccentricity squared

        double g_  [N][N]{};
        double h_  [N][N]{};
        double gSV_[N][N]{};
        double hSV_[N][N]{};
        double baseEpoch_{};
    };

} // namespace sensors

#endif // SAT_PERFORM_IGRF_H