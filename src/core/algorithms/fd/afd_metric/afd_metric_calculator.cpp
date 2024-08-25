#include "afd_metric_calculator.h"

#include <cmath>
#include <iterator>

#include <easylogging++.h>

namespace algos::afd_metric_calculator {

long double AFDMetricCalculator::CalculateG2(config::IndicesType const& lhs_indices,
                                             config::IndicesType const& rhs_indices) const {
    assert(lhs_indices.size() > 0);
    assert(rhs_indices.size() > 0);

    auto num_rows = relation_->GetNumRows();
    assert(num_rows > 0);

    auto num_error_rows = 0.L;

    std::shared_ptr<model::PLI const> lhs_pli = relation_->CalculatePLI(lhs_indices);  // X
    std::shared_ptr<model::PLI const> rhs_pli = relation_->CalculatePLI(rhs_indices);  // Y

    auto const& lhs_clusters = lhs_pli->GetIndex();
    auto pt_shared = rhs_pli->CalculateAndGetProbingTable();
    auto const& pt = *pt_shared.get();
    for (auto const& cluster : lhs_clusters) {
        auto frequencies = model::PLI::CreateFrequencies(cluster, pt);
        auto size = cluster.size();
        if (frequencies.size() != 1 || frequencies.begin()->second != size) num_error_rows += size;
    }

    return num_error_rows / num_rows;
}

std::pair<long double, long double> AFDMetricCalculator::CalculateP1P2(
        size_t num_rows, std::deque<model::PositionListIndex::Cluster>&& lhs_clusters,
        std::deque<model::PositionListIndex::Cluster>&& rhs_clusters) const {
    assert(num_rows > 0);

    auto p1 = 0.L;
    for (auto& y : rhs_clusters) {
        auto size = y.size();
        p1 += size * size;
        std::sort(y.begin(), y.end());
    }
    p1 /= num_rows * num_rows;

    auto p2 = 0.L;
    for (auto& x : lhs_clusters) {
        std::sort(x.begin(), x.end());
        for (auto const& y : rhs_clusters) {
            model::PositionListIndex::Cluster xy;
            std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(xy));

            auto size = (long double)xy.size();
            if (size == 0.L) continue;
            p2 += size * size / x.size();
        }
    }
    p2 /= num_rows;

    return std::make_pair(p1, p2);
}

long double AFDMetricCalculator::CalculateTau(config::IndicesType const& lhs_indices,
                                              config::IndicesType const& rhs_indices) const {
    assert(lhs_indices.size() > 0);
    assert(rhs_indices.size() > 0);

    auto num_rows = relation_->GetNumRows();
    assert(num_rows > 0);

    std::shared_ptr<model::PLI const> rhs_pli = relation_->CalculatePLI(rhs_indices);  // Y
    if (rhs_pli->GetNumCluster() < 2) {
        return 0.L;
    }
    std::shared_ptr<model::PLI const> lhs_pli = relation_->CalculatePLI(lhs_indices);  // X

    auto [p1, p2] = CalculateP1P2(num_rows, lhs_pli->GetAllClusters(), rhs_pli->GetAllClusters());

    return (p2 - p1) / (1 - p1);
}

long double AFDMetricCalculator::CalculateMuPlus(config::IndicesType const& lhs_indices,
                                                 config::IndicesType const& rhs_indices) const {
    assert(lhs_indices.size() > 0);
    assert(rhs_indices.size() > 0);

    auto num_rows = relation_->GetNumRows();
    assert(num_rows > 0);

    std::shared_ptr<model::PLI const> rhs_pli = relation_->CalculatePLI(rhs_indices);  // Y
    if (rhs_pli->GetNumCluster() < 2) {
        return 0.L;
    }
    std::shared_ptr<model::PLI const> lhs_pli = relation_->CalculatePLI(lhs_indices);  // X

    auto lhs_clusters = lhs_pli->GetAllClusters();
    auto x_domain = lhs_clusters.size();
    if (num_rows == x_domain) {
        return 0.L;
    }

    auto [p1, p2] = CalculateP1P2(num_rows, std::move(lhs_clusters), rhs_pli->GetAllClusters());

    long double mu = 1 - (1 - p2) / (1 - p1) * (num_rows - 1) / (num_rows - x_domain);

    return std::max(mu, 0.L);
}

long double AFDMetricCalculator::CalculateFI(config::IndicesType const& lhs_indices,
                                             config::IndicesType const& rhs_indices) const {
    assert(lhs_indices.size() > 0);
    assert(rhs_indices.size() > 0);

    auto num_rows = relation_->GetNumRows();
    assert(num_rows > 0);

    std::shared_ptr<model::PLI const> rhs_pli = relation_->CalculatePLI(rhs_indices);  // Y
    if (rhs_pli->GetNumCluster() < 2) {
        return 0.L;
    }

    auto entropy = rhs_pli->GetEntropy();

    auto rhs_clusters = rhs_pli->GetAllClusters();
    for (auto& y : rhs_clusters) {
        std::sort(y.begin(), y.end());
    }

    auto conditional_entropy = 0.L;
    std::shared_ptr<model::PLI const> lhs_pli = relation_->CalculatePLI(lhs_indices);  // X
    for (auto& x : lhs_pli->GetAllClusters()) {
        std::sort(x.begin(), x.end());
        auto log_x = std::log(x.size());
        for (auto const& y : rhs_clusters) {
            model::PositionListIndex::Cluster xy;
            std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(xy));

            auto size = (long double)xy.size();
            if (size == 0.L) continue;
            conditional_entropy -= size * (std::log(size) - log_x);
        }
    }
    conditional_entropy /= num_rows;

    auto mutual_information = entropy - conditional_entropy;

    return mutual_information / entropy;
}

}  // namespace algos::afd_metric_calculator
