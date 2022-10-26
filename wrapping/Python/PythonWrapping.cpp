#include <PSSolver.hpp>
#include <PSSolver3D.hpp>
#include <PSSolver2D.hpp>
#include <PSSolver1D.hpp>

#include <detail/common.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

namespace py = pybind11;

template <class T> using Holder = std::shared_ptr<T>;

class PyPSSolver : public PSSolver {
public:

  using PSSolver::PSSolver;

};


PYBIND11_MODULE(MODULE_NAME, m) {
  // IMPLEMENT THE WRAPPING HERE
py::class_<PSSolver, PyPSSolver, Holder<PSSolver>>(m, "PSSolver")
      // constructors
    .def(py::init<>())
    .def("calculateBuiltIn", &PSSolver::calculateBuiltIn)
    .def("solveDDQuasiFermi", &PSSolver::solveDDQuasiFermi)
    .def("solveDDClassic", &PSSolver::solveDDClassic)
    .def("setEpsilonGummel", &PSSolver::setEpsilonGummel)
    .def("setPermitivity", &PSSolver::setPermitivity)
    .def("setTemperatur", &PSSolver::setTemperatur)
    .def("setEpsilonNewton", &PSSolver::setEpsilonNewton)
    .def("setElecDensityofStates", &PSSolver::setElecDensityofStates)
    .def("setHolesDensityofStates", &PSSolver::setHolesDensityofStates)
    .def("addSemiconductorSegment", &PSSolver::addSemiconductorSegment)
    .def("addContactSegment", &PSSolver::addContactSegment)
    .def("setCellGrid", &PSSolver::setCellGrid)
    .def("setSimulationGrid", &PSSolver::setSimulationGrid)
    .def("setDonorConcentration", &PSSolver::setDonorConcentration)
    .def("setAcceptorConcentration", &PSSolver::setAcceptorConcentration)
    .def("setNumberGummelIterations", &PSSolver::setNumberGummelIterations)
    .def("setNumberNewtonIterations", &PSSolver::setNumberNewtonIterations)
    .def("addContactPoints", &PSSolver::addContactPoints)
    .def("setValenceBandEnergy", &PSSolver::setValenceBandEnergy)
    .def("setConductionBandEnergy", &PSSolver::setConductionBandEnergy)
    .def("getPotential", &PSSolver::getPotential)
    .def("getElectronConcentration", &PSSolver::getElectronConcentration)
    .def("getHoleConcentration", &PSSolver::getHoleConcentration)
    .def("getChargeDensity", &PSSolver::getChargeDensity)
    .def("setSHR", &PSSolver::setSHR)
    .def("setAuger", &PSSolver::setAuger)
    .def("setSpond", &PSSolver::setSpond)
    .def("setElectronMobility", &PSSolver::setElectronMobility)
    .def("setHoleMobility", &PSSolver::setHoleMobility)
    .def("setSpontanRecomb", &PSSolver::setSpontanRecomb)
    .def("setCn", &PSSolver::setCn)
    .def("setCp", &PSSolver::setCp)
    .def("setTrapEnergy", &PSSolver::setTrapEnergy)
    .def("setElecLifetime", &PSSolver::setElecLifetime)
    .def("setHoleLifetime", &PSSolver::setHoleLifetime)
    .def("applyVoltageSegment", &PSSolver::applyVoltageSegment)
    .def("getEfn", &PSSolver::getEfn)
    .def("getEfp", &PSSolver::getEfp)
    .def("setSteps", &PSSolver::setSteps)
    .def("setChargeDensity", &PSSolver::setChargeDensity)
    .def("solvePoisson", &PSSolver::solvePoisson)
    .def("applyVoltagePoint", &PSSolver::applyVoltagePoint)
    .def("addNeumannSegment", &PSSolver::addNeumannSegment)
    .def("addNeumannPoint", &PSSolver::addNeumannPoint)
    .def("deleteVoltagePoint", &PSSolver::deleteVoltagePoint)
    .def("deleteVoltageSegment", &PSSolver::deleteVoltageSegment); 

py::class_<PSSolver2D, PSSolver, Holder<PSSolver2D>>(m, "PSSolver2D")
    .def(py::init<>())
    .def("import_rectangular_mesh", &PSSolver2D::import_rectangular_mesh)
    .def("test", &PSSolver2D::test);

py::class_<PSSolver1D, PSSolver, Holder<PSSolver1D>>(m, "PSSolver1D")
    .def(py::init<>());

py::class_<PSSolver3D, PSSolver, Holder<PSSolver3D>>(m, "PSSolver3D")
    .def(py::init<>())
    .def("test", &PSSolver3D::test)
    .def("set_polygonpoints", &PSSolver3D::set_polygonpoints)
    .def("set_surfacepolygon", &PSSolver3D::set_surfacepolygon)
    .def("add_simple_boundary", &PSSolver3D::add_simple_boundary);
}


