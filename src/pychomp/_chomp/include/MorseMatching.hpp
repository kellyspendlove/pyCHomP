/// MorseMatching.hpp
/// Shaun Harker
/// 2018-02-23
/// MIT LICENSE

#include "MorseMatching.h"
//#include "EquivariantCubicalMorseMatching.h"
#include "CubicalMorseMatching.h"
#include "CubicalNFMorseMatching.h"
#include "CubicalNNMorseMatching.h"
#include "GenericMorseMatching.h"

inline
std::shared_ptr<MorseMatching>
MorseMatching::compute_matching ( std::shared_ptr<Complex> complex ) {
  if ( std::dynamic_pointer_cast<CubicalComplex>(complex) ) {
    return std::make_shared<CubicalNNMorseMatching>(std::dynamic_pointer_cast<CubicalComplex>(complex));
  } else {
    return std::make_shared<GenericMorseMatching>(complex);
  }
}

inline
std::shared_ptr<MorseMatching>
MorseMatching::compute_matching ( std::shared_ptr<GradedComplex> graded_complex ) {
  if ( std::dynamic_pointer_cast<CubicalComplex>(graded_complex->complex()) ) {
    return std::make_shared<CubicalNNMorseMatching>(graded_complex);
  } else {
    return std::make_shared<GenericMorseMatching>(graded_complex);
  }
}
