#include <gtest/gtest.h>

#include "all_csv_configs.h"
#include "csv_config_util.h"
#include "fd/afd_metric/afd_metric_calculator.h"
#include "indices/type.h"

namespace tests {

struct AFDMetricCalculatorParams {
    CSVConfig csv_config;
    config::IndicesType lhs_indices;
    config::IndicesType rhs_indices;
    long double const tau = 0.L;
    long double const g2 = 0.L;
    long double const fi = 0.L;
    long double const mu_plus = 0.L;

    AFDMetricCalculatorParams(config::IndicesType lhs_indices, config::IndicesType rhs_indices,
                              long double const tau = 0.L, long double g2 = 0.L,
                              long double fi = 0.L, long double mu_plus = 0.L,
                              CSVConfig const& csv_config = kTestFD)
        : csv_config(csv_config),
          lhs_indices(std::move(lhs_indices)),
          rhs_indices(std::move(rhs_indices)),
          tau(tau),
          g2(g2),
          fi(fi),
          mu_plus(mu_plus) {}
};

class TestAFDMetrics : public ::testing::TestWithParam<AFDMetricCalculatorParams> {};

TEST_P(TestAFDMetrics, DefaultTest) {
    auto const& p = GetParam();
    auto input_table = MakeInputTable(p.csv_config);
    std::shared_ptr<ColumnLayoutRelationData> relation =
            ColumnLayoutRelationData::CreateFrom(*input_table, true);
    auto calculator = std::make_unique<algos::afd_metric_calculator::AFDMetricCalculator>(relation);
    EXPECT_DOUBLE_EQ(calculator->CalculateTau(p.lhs_indices, p.rhs_indices), p.tau);
    EXPECT_DOUBLE_EQ(calculator->CalculateG2(p.lhs_indices, p.rhs_indices), p.g2);
    EXPECT_DOUBLE_EQ(calculator->CalculateFI(p.lhs_indices, p.rhs_indices), p.fi);
    EXPECT_DOUBLE_EQ(calculator->CalculateMuPlus(p.lhs_indices, p.rhs_indices), p.mu_plus);
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(
        AFDMetricCalculatorTestSuite, TestAFDMetrics,
        ::testing::Values(
            AFDMetricCalculatorParams({4}, {3}, 78.L/90, 1.L/6, 1 - std::log(4) / std::log(746496), 498.L/630),
            AFDMetricCalculatorParams({3}, {4}, 54.L/114, 5.L/6, std::log(432) / std::log(13824), 252.L/912)
            ));
// clang-format on

}  // namespace tests
