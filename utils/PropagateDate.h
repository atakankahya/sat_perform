//
// Created by ataka on 10.03.2026.
//

#ifndef SAT_PERFORM_PROPAGATEDATE_H
#define SAT_PERFORM_PROPAGATEDATE_H

#include <array>
#include <cmath>
#include <bits/chrono.h>

// Input 0 year 1 month 2 day 3 hour 4 minute 5 second

namespace utils {

    inline std::array<double, 6> propagateDate(const std::array<double, 6>& startDate,
                                               double elapsedSeconds)
    {
        double year = startDate[0];
        double month = startDate[1];
        double day = startDate[2];
        double hour = startDate[3];
        double minute = startDate[4];
        double second = startDate[5];

        double total_second = second + elapsedSeconds;
        second = std::fmod(total_second, 60);

        double total_minutes = minute + std::floor(total_second / 60);
        minute = std::fmod(total_minutes, 60);

        double total_hours = hour + std::floor(total_minutes / 60);
        hour = std::fmod(total_hours, 24);

        double total_days = day + std::floor(total_hours / 24);


        bool is_date_valid = false;

        while (!is_date_valid) {
            int days_in_month= {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

            // check for leap year
            int y = static_cast<int>(year);

            bool is_leap_year = (y%4 == 0 && y%100 != 0) || y%400 == 0;

            if (is_leap_year) {
                days_in_month[1] = 29;
            }

            int current_month_days = days_in_month[static_cast<int>(month) - 1];

            if (total_days > current_month_days) {
                total_days -= current_month_days;
                month = month +1;

                if (month > 12) {
                    month = 1;
                    year += 1;
                }
            } else {
                is_date_valid = true;
            }
        }

        day = total_days;
        return {year, month, day, hour, minute, second};


    }




}
#endif //SAT_PERFORM_PROPAGATEDATE_H