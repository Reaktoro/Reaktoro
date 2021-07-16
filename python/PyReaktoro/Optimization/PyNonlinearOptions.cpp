// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Optimization/NonlinearSolver.hpp>

namespace Reaktoro {

void exportNonlinearOptions(py::module& m)
{
    py::class_<NonlinearOutput, OutputterOptions>(m, "NonlinearOutput")
        .def(py::init<>())
        .def_readwrite("xprefix", &NonlinearOutput::xprefix)
        .def_readwrite("fprefix", &NonlinearOutput::fprefix)
        .def_readwrite("xnames", &NonlinearOutput::xnames)
        .def_readwrite("fnames", &NonlinearOutput::fnames)
        ;

    py::class_<NonlinearOptions>(m, "NonlinearOptions")
        .def(py::init<>())
        .def_readwrite("tolerance", &NonlinearOptions::tolerance)
        .def_readwrite("tolerancex", &NonlinearOptions::tolerancex)
        .def_readwrite("max_iterations", &NonlinearOptions::max_iterations)
        .def_readwrite("tau", &NonlinearOptions::tau)
        .def_readwrite("armijo", &NonlinearOptions::armijo)
        .def_readwrite("linesearch", &NonlinearOptions::linesearch)
        .def_readwrite("output", &NonlinearOptions::output)
        ;
}

} // namespace Reaktoro
