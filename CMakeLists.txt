cmake_minimum_required(VERSION 3.16)
project(sat_perform)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/eigen-5.0.0
)

set(SIM_SOURCES
        spacecraft/SpacecraftState.cpp
        spacecraft/Kinematics.cpp
        sensors/IGRF.cpp
        sensors/Magnetosensor.cpp
        orbit/OrbitPropagator.cpp
        controller/BdotController.cpp
        SatEnv.cpp
)

set(SIM_HEADERS
        SatEnv.h
        config/OrbitConfig.h
        config/SatelliteConfig.h
        config/SimConfig.h
        controller/BdotController.h
        orbit/OrbitPropagator.h
        sensors/IGRF.h
        sensors/Magnetosensor.h
        spacecraft/Kinematics.h
        spacecraft/SpacecraftState.h
        utils/Conversions.h
        utils/GMST.h
        utils/PropagateDate.h
)

add_executable(sat_perform main.cpp ${SIM_SOURCES} ${SIM_HEADERS})
target_include_directories(sat_perform PRIVATE ${PROJECT_INCLUDE_DIRS})

add_executable(debug_magnetic
        spacecraft/SpacecraftState.cpp
        spacecraft/Kinematics.cpp
        sensors/IGRF.cpp
        sensors/Magnetosensor.cpp
        orbit/OrbitPropagator.cpp
        ${SIM_HEADERS}
)
target_include_directories(debug_magnetic PRIVATE ${PROJECT_INCLUDE_DIRS})

include(FetchContent)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.13.6
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(sat_sim bindings.cpp ${SIM_SOURCES} ${SIM_HEADERS})
target_include_directories(sat_sim PRIVATE ${PROJECT_INCLUDE_DIRS})
target_compile_definitions(sat_sim PRIVATE -DNDEBUG)
target_compile_options(sat_sim PRIVATE -O2)