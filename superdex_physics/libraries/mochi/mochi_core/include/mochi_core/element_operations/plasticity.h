#pragma once
#include <mochi_core/utils/batch_types.h>
#include <mochi_core/utils/matrix_utils.h>
#include <cmath>
#include <algorithm>

namespace mochi::fem {
// Local plastic flow keeps the reference metric and simulation topology unchanged.
/// <summary>Transfers only curvature beyond yield into persistent rest curvature, once per accepted step.</summary>
inline bool UpdatePlasticCurvature(
    BatchReal2x2<1> const& referenceMetric,
    BatchReal2x2<1> const& currentCurvatureDifference,
    real yieldCurvature, real flowRate, real stepSeconds,
    BatchReal2x2<1>& plasticCurvature) {
  using V = BatchReal<1>;
  // Measuring with the reference metric avoids mesh-size and rigid-rotation dependent yielding.
  real const determinant = Det(referenceMetric)[0];
  if (!(determinant > 0_r) || !std::isfinite(determinant)) return false;
  auto const elasticDifference = currentCurvatureDifference - plasticCurvature;
  auto const elasticCurvature = Dot(Invert(referenceMetric), elasticDifference);
  real const curvatureMagnitude = std::sqrt(std::max(0_r, Trace(Dot(elasticCurvature, elasticCurvature))[0]));
  if (!std::isfinite(curvatureMagnitude) || curvatureMagnitude <= yieldCurvature) return false;
  // Exponential flow is bounded, dissipative and uses physical elapsed time instead of iteration count.
  real const flowFraction = -std::expm1(-flowRate * stepSeconds);
  real const plasticFraction = flowFraction * (1_r - yieldCurvature / curvatureMagnitude);
  plasticCurvature = plasticCurvature + V{plasticFraction} * elasticDifference;
  return plasticFraction > 0_r;
}
} // namespace mochi::fem
