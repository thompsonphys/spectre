// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"

/// \cond
class ComplexDataVector;
class DataVector;
/// \endcond

namespace GrSelfForce::Tags {

struct MModeRe : db::SimpleTag {
  using type = tnsr::ii<DataVector, 2>;
};

struct MModeIm : db::SimpleTag {
  using type = tnsr::ii<DataVector, 2>;
};

struct Alpha : db::SimpleTag {
  using type = Scalar<ComplexDataVector>;
};

struct Beta : db::SimpleTag {
  using type = Scalar<ComplexDataVector>;
};

struct Gamma : db::SimpleTag {
  using type = tnsr::i<ComplexDataVector, 2>;
};

struct FieldIsRegularized : db::SimpleTag {
  using type = bool;
};

struct SingularField : db::SimpleTag {
  using type = tnsr::ii<ComplexDataVector, 2>;
};

struct BoyerLindquistRadius : db::SimpleTag {
  using type = Scalar<DataVector>;
};

}  // namespace ScalarSelfForce::Tags
