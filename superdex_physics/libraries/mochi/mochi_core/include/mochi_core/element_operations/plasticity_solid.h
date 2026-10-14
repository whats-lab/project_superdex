#pragma once
#include <mochi_core/utils/batch_types.h>
#include <mochi_core/utils/decomposition_utils.h>
#include <mochi_core/utils/matrix_utils.h>
#include <mochi_core/utils/span.h>
#include <cmath>
#include <algorithm>

namespace mochi::fem {
// Plastic reference data stays fixed while the nonlinear solver evaluates candidate positions.
/// <summary>Maps energy, stress and tangent through Fe = F times inverse plastic deformation.</summary>
template<int kBatchSize, class Response>
auto PlasticConstitutiveResponse(Response const& response, Span<NdArray<real,3,3> const> inversePlastic) {
  return [&response, inversePlastic](NdArray<int,kBatchSize> const& indices,
      BatchReal3x3<kBatchSize> const& deformation,
      BatchDouble<kBatchSize>* energy, BatchReal3x3<kBatchSize>* stress,
      NdArray<BatchReal3x3<kBatchSize>,3,3>* tangent, bool projectPsd) {
    using V = BatchReal<kBatchSize>;
    // Identity/OFF uses the original material callback without transformed arithmetic.
    if(inversePlastic.empty()) { response(indices,deformation,energy,stress,tangent,projectPsd); return; }
    BatchReal3x3<kBatchSize> inverse;
    for(int row=0;row<3;++row) for(int column=0;column<3;++column) {
      alignas(alignof(V)) real values[V::kSize]{};
      for(int lane=0;lane<kBatchSize;++lane) values[lane]=inversePlastic[indices[lane]][row][column];
      inverse[row][column]=Load<V>(values);
    }
    auto elastic=Dot(deformation,inverse);
    response(indices,elastic,energy,stress,tangent,projectPsd);
    // Chain-rule transforms preserve the gradient and symmetric material tangent in original coordinates.
    if(stress) *stress=Dot(*stress,Transpose(inverse));
    if(tangent) {
      auto original=*tangent;
      for(int row=0;row<3;++row) for(int column=0;column<3;++column)
        for(int otherRow=0;otherRow<3;++otherRow) for(int otherColumn=0;otherColumn<3;++otherColumn) {
          V value{0_r};
          for(int inner=0;inner<3;++inner) for(int otherInner=0;otherInner<3;++otherInner)
            value+=original[row][inner][otherRow][otherInner]*inverse[column][inner]*inverse[otherColumn][otherInner];
          (*tangent)[row][column][otherRow][otherColumn]=value;
        }
    }
  };
}

/// <summary>Returns a rotation-invariant plastic reference candidate from yield-exceeding principal log strains.</summary>
inline bool PlasticDeformationCandidate(BatchReal3x3<1> const& deformation,
    BatchReal3x3<1> const& inversePlastic, real yieldStrain, real flowRate, real seconds,
    BatchReal3x3<1>& candidate) {
  using V=BatchReal<1>;
  auto const elastic=Dot(deformation,inversePlastic);
  // An inverted or collapsed element cannot define a valid new material reference.
  if(!(Det(elastic)[0]>1e-12_r)) return false;
  BatchReal3x3<1> left,rightTranspose;
  BatchReal3<1> stretch;
  BatchedRotationVariantSvd3x3<1>(elastic,left,stretch,rightTranspose);
  real logarithms[3], squaredMagnitude=0_r;
  for(int axis=0;axis<3;++axis) {
    real const value=stretch[axis][0];
    if(!(value>1e-6_r)||!std::isfinite(value)) return false;
    logarithms[axis]=std::log(value); squaredMagnitude+=logarithms[axis]*logarithms[axis];
  }
  real const magnitude=std::sqrt(squaredMagnitude);
  if(!std::isfinite(magnitude)||magnitude<=yieldStrain) return false;
  real const fraction=-std::expm1(-flowRate*seconds)*(1_r-yieldStrain/magnitude);
  auto diagonal=Eye<3,V>();
  for(int axis=0;axis<3;++axis) diagonal[axis][axis]=V{std::exp(-fraction*logarithms[axis])};
  candidate=Dot(inversePlastic,Dot(Transpose(rightTranspose),Dot(diagonal,rightTranspose)));
  return Det(candidate)[0]>1e-12_r && IsFinite(candidate);
}
} // namespace mochi::fem
