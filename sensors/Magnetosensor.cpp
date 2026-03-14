//
// Created by ataka on 2.03.2026.
//

#include "Magnetosensor.h"
#include "../utils/Conversions.h"
#include <iostream>


namespace sensors {
    MagnetoSensor::MagnetoSensor(const IGRF &igrf)
        : igrf_(igrf)
    {
    }

    Eigen::Vector3d MagnetoSensor::measure(const orbit::LLA &lla, double gmst_rad, double decimal_year, const Eigen::Matrix3d &R_eci2body) const {

        // b in ned from the igrf model
        Eigen::Vector3d B_ned = igrf_.computeNED(lla.lat_deg,lla.lon_deg,lla.alt_km,decimal_year);
        //std::cout << B_ned(0)<< "\n" << B_ned(1) << "\n" << B_ned(2) << std::endl;
        // rotation chain ned to ecef to eci to body
        Eigen::Matrix3d R_ned2ecef = utils::nedToECEFDCM(lla.lat_deg,lla.lon_deg);
        Eigen::Matrix3d R_ecef2eci = utils::ecefToEciDCM(gmst_rad);

        // R_eci2body must map vectors expressed in ECI into BODY.
        //Eigen::Vector3d a = R_eci2body * R_ecef2eci * R_ned2ecef * B_ned;
        //std::cout << "Bbody" << a(0) << "  " <<a(1) << "  " <<a(2) << std::endl;
        return R_eci2body * R_ecef2eci * R_ned2ecef * B_ned;



    }
}
