cmake_minimum_required(VERSION 4.1)
project(sat_perform)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(eigen-5.0.0)
add_executable(sat_perform main.cpp
        spacecraft/SpacecraftState.cpp
        spacecraft/SpacecraftState.h
        sensors/Magnetosensor.cpp
        sensors/Magnetosensor.h
        spacecraft/Kinematics.cpp
        spacecraft/Kinematics.h
        utils/Conversions.h
        config/SatelliteConfig.h
        config/SimConfig.h
        config/OrbitConfig.h
        orbit/OrbitPropagator.cpp
        orbit/OrbitPropagator.h
        utils/PropagateDate.h
        sensors/IGRF.cpp
        sensors/IGRF.h)
