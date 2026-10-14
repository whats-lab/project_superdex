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
  std::vector<NdArray<real, 2, 2>> curvature;
  MOCHI_STRUCT_BEGIN(mochi::CShellPlasticity)
  MOCHI_ATTRIBUTE(CaptureState)
  MOCHI_FIELD(enabled)
  MOCHI_FIELD(yieldCurvature)
  MOCHI_FIELD(flowRate)
  MOCHI_FIELD(curvature)
  MOCHI_STRUCT_END()

};
} // namespace mochi
