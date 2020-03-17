/// CubicalMorseMatching.h
/// Shaun Harker
/// 2018-02-16
/// MIT LICENSE

#pragma once

#include <memory>
#include <unordered_set>
#include <vector>

#include "Integer.h"
#include "Chain.h"
#include "Complex.h"
#include "GradedComplex.h"
#include "MorseMatching.h"

class EquivariantCubicalMorseMatching : public MorseMatching {
public:
  /// CubicalMorseMatching
  EquivariantCubicalMorseMatching ( std::shared_ptr<CubicalComplex> complex_ptr ) : complex_(complex_ptr) {
    type_size_ = complex_ -> type_size();
    graded_complex_.reset(new GradedComplex(complex_, [](Integer i){return 0;}));
  }

  /// CubicalMorseMatching
  EquivariantCubicalMorseMatching ( std::shared_ptr<GradedComplex> graded_complex_ptr ) : graded_complex_(graded_complex_ptr) {
    complex_ = std::dynamic_pointer_cast<CubicalComplex>(graded_complex_->complex());
    if ( not complex_ ) {
      throw std::invalid_argument("CubicalMorseMatching must be constructed with a Cubical Complex");
    }
    type_size_ = complex_ -> type_size();
    Integer D = complex_ -> dimension();
    Integer idx = 0;
    begin_.resize(D+2);
    for ( Integer d = 0; d <= D; ++ d) {
      begin_[d] = idx;
      for ( auto v : (*complex_)(d) ) { // TODO: skip fringe cells
        if ( ! complex_ -> rightfringe(v) ) {
          if ( mate(v) == v ) { 
            reindex_.push_back({v,idx});
            ++idx;
          }
        }
      }
    }
    begin_[D+1] = idx;
  }

  /// critical_cells
  std::pair<BeginType const&,ReindexType const&>
  critical_cells ( void ) const {
    return {begin_,reindex_};
  }
  std::vector<Integer> 
  permute_tupled_barys ( std::vector<std::tuple<Integer,Integer>> tuples, std::vector<Integer> permutation ) const {
    //Given vector of tuples of barycenter coordinates and permutation
    //Apply permutation to vector of tuples
    //Hand back untupled barycenter coordinates
    std::vector <Integer> permuted_barys;
    for (Integer i = 0; i < tuples . size(); i++ ){
      permuted_barys . push_back ( std::get<0>(tuples[permutation[i]]) );
      permuted_barys . push_back ( std::get<1>(tuples[permutation[i]]) );
    }
    return permuted_barys;
  }
  std::vector<Integer> 
  permute_barys ( std::vector<Integer> barys, std::vector<Integer> &permuted_barys ) const { 

    //for ( auto y : barys ) std::cout << y << ", "; std::cout << "]\n";
    
    std::vector <std::tuple<Integer, Integer, Integer>> tuples_to_sort;
    for (Integer i = 0; i < barys.size(); i+= 2 ) {
      tuples_to_sort . push_back ( std::make_tuple(barys[i],barys[i+1], i/2));
    }
    std::sort(tuples_to_sort.begin(),tuples_to_sort.end());

    //for ( auto y : tuples_to_sort ) std::cout << std::get<0>(y) << " " << std::get<1>(y) << " ,"; std::cout << "\n";

    //std::vector <std::tuple<Integer, Integer>> new_tuples;
    std::vector <Integer> permutation;
    for (Integer i = 0; i < tuples_to_sort.size(); i++ ) {
      permuted_barys . push_back ( std::get<0>(tuples_to_sort[i]) );
      permuted_barys . push_back ( std::get<1>(tuples_to_sort[i]) );
      //new_tuples . push_back ( std::make_tuple(std::get<0>(tuples_to_sort[i]), std::get<1>(tuples_to_sort[i])) );
      permutation . push_back ( std::get<2>(tuples_to_sort[i]));
    }
    Integer count = 0;
    std::vector<Integer> unpermute;
    unpermute . resize (permutation.size() );
    for ( auto y : permutation ) {
      unpermute [y] = count;
      ++count;
    }
    //for ( auto y : permuted_barys ) std::cout << y << ", "; std::cout << "]\n";
    //for ( auto y : permuted_barys ) std::cout << " ---- "; std::cout << "]\n";

    return unpermute;
  }
  Integer
  unpermute_index ( Integer x, std::vector<Integer> unperm) const {
    // Create tuples of coords
    std::vector<Integer> barys = complex_ -> barycenter(x);
    std::vector <std::tuple<Integer, Integer>> tuples;
    for (Integer i = 0; i < barys.size(); i+= 2 ) {
      tuples . push_back ( std::make_tuple(barys[i],barys[i+1]) );
    }
    //call permute tupled coords with unpermute, gives back flattened
    std::vector<Integer> unpermuted_barys = permute_tupled_barys ( tuples, unperm );
    return complex_ -> IB(unpermuted_barys);
  }
  // std::vector<Integer> 
  // permute_barys ( std::vector<Integer> barys, std::vector<Integer> &upermute ) const { 

  //   for ( auto y : barys ) std::cout << y << ", "; std::cout << "]\n";
    
  //   std::vector <std::tuple<Integer, Integer, Integer>> tupled_barys;
  //   for (Integer i = 0; i < barys.size(); i+= 2 ) {
  //     tupled_barys . push_back ( std::make_tuple(barys[i],barys[i+1], i/2));
  //   }
  //   std::sort(tupled_barys.begin(),tupled_barys.end());

  //   std::vector <std::tuple<Integer, Integer>> new_tuples;
  //   std::vector <Integer> permutation;
  //   for (Integer i = 0; i < tupled_barys.size(); i++ ) {
  //     new_tuples . push_back ( std::make_tuple(std::get<0>(tupled_barys[i]), std::get<1>(tupled_barys[i])) );
  //     permutation . push_back ( std::get<2>(tupled_barys[i]));
  //   }
  //   std::vector <Integer> permuted_barys;
  //   for (Integer i = 0; i < new_tuples . size(); i++ ){
  //     permuted_barys . push_back ( std::get<0>(new_tuples[permutation[i]]) );
  //     permuted_barys . push_back ( std::get<1>(new_tuples[permutation[i]]) );
  //   }

  //   for ( auto y : permuted_barys ) std::cout << y << ", "; std::cout << "]\n";
  //   for ( auto y : permuted_barys ) std::cout << " ---- "; std::cout << "]\n";

  //   return permuted_barys;
  // }
  /// mate
  Integer
  mate ( Integer x ) const {
    //Given index, hand to permute_barys to get permuted barycenter coordinates, and the permutation
    //Give permuted barycenter coordinates to IB to get the permuted index
    //Apply mate to index
    //Unpermute mate
    std::vector<Integer> permuted_barys;
    std::vector<Integer> unpermute = permute_barys(complex_ -> barycenter(x), permuted_barys);
    Integer permuted_index = complex_ -> IB(permuted_barys);
    Integer permuted_mate = mate_(permuted_index, complex_ -> dimension() );
    Integer unpermuted_mate = unpermute_index ( permuted_mate, unpermute );
    return unpermuted_mate;
    //Integer unpermuted_mate = (perm)
    //return permute_index ( mate_(permute_index(x), complex_ -> dimension() ) );
    
    //return mate_(x, complex_ -> dimension());
  }

  /// priority
  Integer
  priority ( Integer x ) const { 
    return type_size_ - x % type_size_;
  }

private:
  uint64_t type_size_;
  std::shared_ptr<GradedComplex> graded_complex_;
  std::shared_ptr<CubicalComplex> complex_;
  BeginType begin_;
  ReindexType reindex_;

  // def mate(cell, D):
  // for d in range(0, D):
  //   if cell has extent in dimension d:
  //     left = leftboundary(cell, d)
  //     if value(left) == value(cell):
  //       if left == mate(left, d):
  //         return left
  //   else:
  //     right = rightcoboundary(cell, d)
  //     if value(right) == value(cell):
  //       if right == mate(right, d):
  //         return right
  //   return cell 
  // Note: should the complicated formulas (which are also found in CubicalComplex.h not be repeated here?
  // Note: the reason for the "fringe" check preventing mating is that otherwise it is possible to 
  //       end up with a cycle 
  // TODO: Furnish a proof of correctness and complexity that this cannot produce cycles.
  Integer mate_ ( Integer cell, Integer D ) const {
    //bool fringe = complex_ -> rightfringe(cell);
    if ( complex_ -> rightfringe(cell) ) return cell; // MAYBE
    //Integer mincoords = complex_ -> mincoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
    //Integer maxcoords = complex_ -> maxcoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
    Integer shape = complex_ -> cell_shape(cell);
    Integer position = cell % complex_ -> type_size();
    if ( position == complex_ -> type_size() - 1 ) return cell; // Break cycles
    for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= 1L  ) {
      // If on right fringe for last dimension, prevent mating with left boundary
      if ( (d == D-1) && (position + complex_ -> PV()[d] >= complex_ -> type_size()) ) break;
      //if ( fringe && (mincoords & bit) ) continue; // Todo -- is this the best
      //if ( bit & maxcoords ) continue; // Don't connect fringe to acyclic part
      Integer type_offset = complex_ -> type_size() * complex_ -> TS() [ shape ^ bit ];
      Integer proposed_mate = position + type_offset;
      if ( graded_complex_ -> value(proposed_mate) == graded_complex_ -> value(cell) && proposed_mate == mate_(proposed_mate, d) ) { 
        return proposed_mate;
      }
    }
    return cell;
  }
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
EquivariantCubicalMorseMatchingBinding(py::module &m) {
  py::class_<EquivariantCubicalMorseMatching, std::shared_ptr<EquivariantCubicalMorseMatching>>(m, "EquivariantCubicalMorseMatching")
    .def(py::init<std::shared_ptr<CubicalComplex>>())
    .def(py::init<std::shared_ptr<GradedComplex>>())    
    .def("mate", &EquivariantCubicalMorseMatching::mate)
    .def("priority", &EquivariantCubicalMorseMatching::priority);
}


// inline void
// CubicalMorseMatchingBinding(py::module &m) {
//   py::class_<CubicalMorseMatching, std::shared_ptr<CubicalMorseMatching>>(m, "CubicalMorseMatching")
//     .def(py::init<std::shared_ptr<CubicalComplex>>())
//     .def(py::init<std::shared_ptr<GradedComplex>>())    
//     .def("mate", &CubicalMorseMatching::mate)
//     .def("priority", &CubicalMorseMatching::priority);
// }