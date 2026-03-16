// Created by ataka on 14.03.2026.
//

#include <pybind11/stl.h>
#include "SatEnv.h"

namespace py = pybind11;

PYBIND11_MODULE(sat_sim, m) {
    m.doc() = "UPMSat-2 magnetic ADCS simulation environment";

    py::class_<SatEnv>(m, "SatEnv")
        .def(py::init<const std::string&>(), py::arg("igrf_datadir"))
        .def("reset", &SatEnv::reset)
        .def("step", &SatEnv::step)
        .def_readonly_static("obs_dim", &SatEnv::OBS_DIM)
        .def_readonly_static("n_actions", &SatEnv::N_ACTIONS);
}
