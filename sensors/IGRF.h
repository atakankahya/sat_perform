//
// Created by ataka on 10.03.2026.
//

#ifndef SAT_PERFORM_IGRF_H
#define SAT_PERFORM_IGRF_H


#include "Eigen/Dense"


namespace sensors {
    class IGRF {
        public:
        explicit IGRF(const std::string& datadir);
        [[nodiscard]] Eigen::Vector3d computeNED(double lat_deg,
                                                double lon_deg,
                                                double alt_km,
                                                double decimal_year) const;

        private:
        static constexpr int NMAX = 13;
        static constexpr int N    = 14;

        static constexpr double RE = 6371.1;
        static constexpr double WGS84_a = 6378.137;
        static constexpr double WGS84_e2 = 0.0069437999014;

        double g_ [N][N]{};
        double h_ [N][N]{};
        double gSV_ [N][N]{};
        double hSV_ [N][N]{};
        double baseEpoch_ {};

    };


}

#endif //SAT_PERFORM_IGRF_H