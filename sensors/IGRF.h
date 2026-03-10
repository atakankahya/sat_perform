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


    };


}

#endif //SAT_PERFORM_IGRF_H