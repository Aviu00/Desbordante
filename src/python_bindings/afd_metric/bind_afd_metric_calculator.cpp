#include "bind_afd_metric_calculator.h"

#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/fd/afd_metric/afd_metric_calculator.h"
#include "descriptions.h"
#include "py_util/py_to_any.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {
using namespace pybind11::literals;

void BindAfdMetricCalculation(py::module_& main_module) {
    using namespace algos;
    using namespace algos::afd_metric_calculator;
    using namespace config::descriptions;

    auto desc_metric_str = [](char const* metric) {
        return std::string("Calculates ") + metric + " metric on specified indices.\n";
    };
    auto indices_input_str = std::string("Inputs:\n\tlhs_indices: ") + kDLhsIndices +
                             "\n\trhs_indices: " + kDRhsIndices;

    auto afd_metric_module = main_module.def_submodule("afd_metric_calculation");
    py::class_<AFDMetricCalculator>(afd_metric_module, "AFDMetricCalculator")
            .def(py::init([](py::handle table, bool is_null_eq_null) {
                     return std::make_unique<AFDMetricCalculator>(PyToInputTable(table),
                                                                  is_null_eq_null);
                 }),
                 "table"_a, "is_null_eq_null"_a = true,
                 (std::string("Inputs:\n\ttable: ") + kDTable +
                  "\n\tis_null_eq_null: " + kDEqualNulls)
                         .c_str())
            .def("calculate_g2", &AFDMetricCalculator::CalculateG2, "lhs_indices"_a,
                 "rhs_indices"_a, (desc_metric_str("G2") + indices_input_str).c_str())
            .def("calculate_tau", &AFDMetricCalculator::CalculateTau, "lhs_indices"_a,
                 "rhs_indices"_a, (desc_metric_str("τ") + indices_input_str).c_str())
            .def("calculate_mu_plus", &AFDMetricCalculator::CalculateMuPlus, "lhs_indices"_a,
                 "rhs_indices"_a, (desc_metric_str("μ+") + indices_input_str).c_str())
            .def("calculate_fi", &AFDMetricCalculator::CalculateFI, "lhs_indices"_a,
                 "rhs_indices"_a, (desc_metric_str("FI") + indices_input_str).c_str());
}
}  // namespace python_bindings
