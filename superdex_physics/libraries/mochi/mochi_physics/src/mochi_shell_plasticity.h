#pragma once
#include "mochi_common_components.h"
#include <mochi_core/utils/nd_array.h>
#include <vector>

namespace mochi {
// Each actor owns its optional history; default-OFF actors allocate no curvature array.
struct CShellPlasticity : NoCopy {
  bool enabled = false;
  real yieldCurvature = 8_r;
  real flowRate = 30_r;
  bool useRegionOverrides = false;
  bool regionsOnly = false;
  std::vector<real> elementYield;
  std::vector<real> elementFlowRate;
  std::vector<NdArray<real, 2, 2>> curvature;
  std::vector<NdArray<real, 2, 2>> edgeCurvature;
  MOCHI_STRUCT_BEGIN(mochi::CShellPlasticity)
  MOCHI_ATTRIBUTE(CaptureState)
  MOCHI_FIELD(enabled)
  MOCHI_FIELD(yieldCurvature)
  MOCHI_FIELD(flowRate)
  MOCHI_FIELD(useRegionOverrides)
  MOCHI_FIELD(regionsOnly)
  MOCHI_FIELD(elementYield)
  MOCHI_FIELD(elementFlowRate)
  MOCHI_FIELD(curvature)
  MOCHI_FIELD(edgeCurvature)
  MOCHI_STRUCT_END()

};
// Thickness changes elastic membrane and bending stiffness even while plasticity is disabled.
struct CShellThickness : NoCopy {
  std::vector<real> ratios;
  MOCHI_STRUCT_BEGIN(mochi::CShellThickness)
  MOCHI_ATTRIBUTE(CaptureState)
  MOCHI_FIELD(ratios)
  MOCHI_STRUCT_END()
};
} // namespace mochi
