#include "rfc_1149_5.hpp"

const std::array<double, 3> pass::rfc_1149_5_search::RANDOM_VECTOR {4, 4, 4};

pass::rfc_1149_5_search::rfc_1149_5_search() noexcept {
  this->optimisation_function = [this](const auto& problem, const auto& _) {
    mant::optimise_result<double, 3> result;
    auto&& start_time = std::chrono::steady_clock::now();

    result.parameter = RANDOM_VECTOR;
    result.objective_value = problem.objective_function(RANDOM_VECTOR);
    ++result.evaluations;
    result.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - start_time);
    return result;
  };
}
