//
// Created by ataka on 2.03.2026.
//

#ifndef SAT_PERFORM_MAGNETOSENSOR_H
#define SAT_PERFORM_MAGNETOSENSOR_H
#include "IGRF.h"
#include "../orbit/OrbitPropagator.h"

namespace sensors {
    class MagnetoSensor {
        public:
        explicit MagnetoSensor(const IGRF& igrf);

        [[nodiscard]] Eigen::Vector3d measure(

            const orbit::LLA& lla,
            double gmst_rad,
            double decimal_year,
            const Eigen::Matrix3d& R_eci2body) const;
    private:
        const IGRF& igrf_;
    };

}



#endif //SAT_PERFORM_MAGNETOSENSOR_H
