#pragma once

#include <memory>

#include "indices/type.h"
#include "table/column_layout_relation_data.h"
#include "tabular_data/input_table_type.h"

namespace algos::afd_metric_calculator {

class AFDMetricCalculator {
private:
    std::shared_ptr<ColumnLayoutRelationData> relation_;

    std::pair<long double, long double> CalculateP1P2(
            size_t num_rows, std::deque<model::PositionListIndex::Cluster>&& lhs_clusters,
            std::deque<model::PositionListIndex::Cluster>&& rhs_clusters) const;

public:
    long double CalculateG2(config::IndicesType const& lhs_indices,
                            config::IndicesType const& rhs_indices) const;

    long double CalculateTau(config::IndicesType const& lhs_indices,
                             config::IndicesType const& rhs_indices) const;

    long double CalculateMuPlus(config::IndicesType const& lhs_indices,
                                config::IndicesType const& rhs_indices) const;

    long double CalculateFI(config::IndicesType const& lhs_indices,
                            config::IndicesType const& rhs_indices) const;

    AFDMetricCalculator(std::shared_ptr<ColumnLayoutRelationData> relation)
        : relation_(std::move(relation)) {};

    AFDMetricCalculator(config::InputTable input_table, bool is_null_eq_null = true)
        : relation_(ColumnLayoutRelationData::CreateFrom(*input_table, is_null_eq_null)) {};
};

}  // namespace algos::afd_metric_calculator
