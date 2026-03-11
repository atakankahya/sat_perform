//
// Created by ataka on 2.03.2026.
//

#include "Magnetosensor.h"
#include "../utils/Conversions.h"



namespace sensors {
    MagnetoSensor::MagnetoSensor(const IGRF &igrf)
        : igrf_(igrf)
    {
    }

    Eigen::Vector3d MagnetoSensor::measure(const orbit::LLA &lla, double gmst_rad, double decimal_year, const Eigen::Matrix3d &R_eci2body) const {

        // b in ned from the igrf model
        Eigen::Vector3d B_ned = igrf_.computeNED(lla.lat_deg,lla.lon_deg,lla.alt_km,decimal_year);

        // rotation chain ned to ecef to eci to body
        Eigen::Matrix3d R_ned2ecef = utils::nedToECEFDCM(lla.lat_deg,lla.lon_deg);
        Eigen::Matrix3d R_ecef2eci = utils::ecefToEciDCM(gmst_rad);


        // b in body frame
        return R_eci2body * R_ned2ecef * B_ned;

    }
}
