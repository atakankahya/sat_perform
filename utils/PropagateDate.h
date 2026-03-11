//
// Created by ataka on 10.03.2026.
//

#ifndef SAT_PERFORM_PROPAGATEDATE_H
#define SAT_PERFORM_PROPAGATEDATE_H

#include <array>
#include <cmath>
#include <bits/chrono.h>
#include <ctime>

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
            int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

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

    struct Date {
        int year;
        int month;
        int day;
        int hour;
        int minute;
        double second;
    };

    inline double dateToDecimalYear(const Date& d) {
        int dpm[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        bool leap = (d.year % 400 == 0) || (d.year % 4 == 0 && d.year % 100 != 0);
        if (leap) dpm[1] = 29;

        int days_in_year = leap ? 366 : 365;

        int doy = 0;
        for (int i = 0; i < d.month - 1; ++i)
            doy += dpm[i];
        doy += d.day - 1;

        double frac_day = (d.hour + 3600.0 + d.minute * 60.0 + d.second) / 86400.0;
        return static_cast<double>(d.year) + (doy + frac_day) / static_cast<double>(days_in_year);
    }

    inline double advanceDecimalYear(double base_decimal_year, double elapsed_seconds) {
        return base_decimal_year + elapsed_seconds / (365.25 * 86400.0);
    }

    double getDecimalYear(int year, int month, int day) {
        std::tm t = {};
        t.tm_year = year - 1900;
        t.tm_mon = month - 1;
        t.tm_mday = day;

        std::mktime(&t);

        int yday = t.tm_yday + 1;
        bool leap = (year % 4== 0 && year % 100 != 0) || year % 400 == 0;
        int days_in_year = leap ? 366 : 365;

        return static_cast<double>(year) + static_cast<double>(yday) / static_cast<double>(days_in_year);
    }




}
#endif //SAT_PERFORM_PROPAGATEDATE_H