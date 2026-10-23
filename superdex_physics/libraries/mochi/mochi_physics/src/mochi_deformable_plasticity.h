#pragma once
#include "mochi_common_components.h"
#include <mochi_core/utils/nd_array.h>
#include <vector>

namespace mochi {
// Empty components keep capture/restore entity sets stable across runtime feature switches.
struct CSoftPlasticity : NoCopy {
  bool enabled=false;
  real yieldStrain=0.1_r;
  real flowRate=30_r;
  bool useRegionOverrides=false;
  std::vector<real> elementYield;
  std::vector<real> elementFlowRate;
  std::vector<NdArray<real,3,3>> inverseDeformation;
  MOCHI_STRUCT_BEGIN(mochi::CSoftPlasticity)
  MOCHI_ATTRIBUTE(CaptureState)
  MOCHI_FIELD(enabled)
  MOCHI_FIELD(yieldStrain)
  MOCHI_FIELD(flowRate)
  MOCHI_FIELD(useRegionOverrides)
  MOCHI_FIELD(elementYield)
  MOCHI_FIELD(elementFlowRate)
  MOCHI_FIELD(inverseDeformation)
  MOCHI_STRUCT_END()
};
struct CRodPlasticity : NoCopy {
  bool enabled=false;
  real yieldStrain=0.1_r;
  real yieldCurvature=8_r;
  real yieldTwist=8_r;
  real flowRate=30_r;
  std::vector<NdArray<real,4>> strains;
  MOCHI_STRUCT_BEGIN(mochi::CRodPlasticity)
  MOCHI_ATTRIBUTE(CaptureState)
  MOCHI_FIELD(enabled)
  MOCHI_FIELD(yieldStrain)
  MOCHI_FIELD(yieldCurvature)
  MOCHI_FIELD(yieldTwist)
  MOCHI_FIELD(flowRate)
  MOCHI_FIELD(strains)
  MOCHI_STRUCT_END()
};
} // namespace mochi
