//
// GMST calculation matching MATLAB implementation
//

#ifndef SAT_PERFORM_GMST_H
#define SAT_PERFORM_GMST_H

#include <cmath>
#include <array>

namespace utils {
    // Calculate Julian Date from date array
    inline double computeJulianDate(const std::array<double, 6>& date) {
        double Y = date[0];
        double Mo = date[1];
        double D = date[2];
        double h = date[3];
        double mi = date[4];
        double s = date[5];

        // adjust for January/February
        if (Mo <= 2) {
            Y = Y - 1;
            Mo = Mo + 12;
        }

        double A = std::floor(Y / 100.0);
        double B = 2.0 - A + std::floor(A / 4.0);

        double JD = std::floor(365.25 * (Y + 4716.0)) +
                    std::floor(30.6001 * (Mo + 1.0)) +
                    D + B - 1524.5 +
                    (h + mi/60.0 + s/3600.0) / 24.0;

        return JD;
    }

    // Calculate GMST in radians from date and elapsed time
    inline double computeGMST(const std::array<double, 6>& start_date, double elapsed_seconds) {
        // Calculate Julian Date at current time
        std::array<double, 6> current_date = start_date;

        // Add elapsed seconds to the date
        double total_seconds = current_date[5] + elapsed_seconds;
        double add_minutes = std::floor(total_seconds / 60.0);
        current_date[5] = std::fmod(total_seconds, 60.0);

        double total_minutes = current_date[4] + add_minutes;
        double add_hours = std::floor(total_minutes / 60.0);
        current_date[4] = std::fmod(total_minutes, 60.0);

        double total_hours = current_date[3] + add_hours;
        double add_days = std::floor(total_hours / 24.0);
        current_date[3] = std::fmod(total_hours, 24.0);

        current_date[2] += add_days;
        // Simplified - doesn't handle month/year rollover, but good enough for short sims

        double JD = computeJulianDate(current_date);

        // Calculate T (centuries from J2000)
        double T = (JD - 2451545.0) / 36525.0;

        // GMST in degrees - IAU formula
        double theta = 280.46061837 +
                      360.98564736629 * (JD - 2451545.0) +
                      0.000387933 * T * T -
                      (T * T * T) / 38710000.0;

        // Normalize to [0, 360)
        theta = std::fmod(theta, 360.0);
        if (theta < 0.0) theta += 360.0;

        // Convert to radians
        return theta * M_PI / 180.0;
    }
}

#endif //SAT_PERFORM_GMST_H