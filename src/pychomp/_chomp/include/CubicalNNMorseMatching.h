/// CubicalNNMorseMatching.h
/// Shaun Harker & Kelly Spendlove
/// 2018-02-16
/// MIT LICENSE

/*
This is an effort to not iterate over fringe cells
*/

#pragma once

#include <memory>
#include <unordered_set>
#include <vector>
#include <queue>

#include "Integer.h"
#include "Chain.h"
#include "Complex.h"
#include "CubicalComplex.h"
#include "GradedComplex.h"
#include "MorseMatching.h"

class CubicalNNMorseMatching : public MorseMatching {
public:
  /// CubicalNNMorseMatching
  CubicalNNMorseMatching ( std::shared_ptr<CubicalComplex> complex_ptr ) : complex_(complex_ptr) {
    type_size_ = complex_ -> type_size();
    graded_complex_.reset(new GradedComplex(complex_, [](Integer i){return 0;}));
  }

  /// CubicalNNMorseMatching
  CubicalNNMorseMatching ( std::shared_ptr<GradedComplex> graded_complex_ptr ) : graded_complex_(graded_complex_ptr) {
    complex_ = std::dynamic_pointer_cast<CubicalComplex>(graded_complex_->complex());
    if ( not complex_ ) {
      throw std::invalid_argument("CubicalNNMorseMatching must be constructed with a Cubical Complex");
    }
    type_size_ = complex_ -> type_size();
    

    Integer D = complex_ -> dimension(); //dimensions for Morse complex
    //Optimizations for Configuration Space (of squares in a rectangle)
    // Integer n = complex_ -> dimension() / 2;
    // Integer p = complex_ -> boxes () [ 0 ];
    // Integer q = complex_ -> boxes () [ n ];
    // D = std::min ( n, std::min( (p*q)/3, p*q-n ));//adjust dimension using bounds from paper
    //End optimizations for squares
    //Optimizations for Configuration Space (of cubes in a box)
    // Integer n = complex_ -> dimension() / 3;
    // Integer p = complex_ -> boxes () [ 0 ];
    // Integer q = complex_ -> boxes () [ n ];
    // Integer w = complex_ -> boxes () [ 2*n ];
    // Integer D = std::min ( 2*n, std::min((p*q*w*3) / 8, p*q*w-n) );
    //End optimizations for cubes

    //auto & ST = complex_ -> ST();
    begin_.resize(D+2);


    std::vector<std::queue<Integer>> crit_cells (D+1);
    for ( auto v : (*complex_)(0) ) { //for each position (vertex)

      if ( v % (complex_ -> PV()[complex_ -> dimension()] / 10) == 0) std::cout << "progress: " << v << "\n";
      //config space optimization: any cell incident to non-configuration doesn't belong to config space
      //if ( graded_complex_ -> value(v) != 0 ) continue; //config optimization

      Integer maxc = complex_ -> maxcoords ( v );
      Integer H = complex_ -> dimension() - popcount_ ( maxc ); //H is size of hypercube to iterate over
      std::vector<bool> visited ( (1L<<H) );

      //Braids optimizer
      //TODO: Furnish argument why this works if maxc > 0
      //if ( maxc == 0 ) {//Avoid right fringe?
      // Integer top_cell = v + complex_ -> type_size() * complex_ -> TS() [ embed((1<<H)-1, maxc) ];
      // // // //if (H != 0 && graded_complex_ -> value (v) == graded_complex_ -> value ( top_cell) ) continue; //MAYBE
      // auto top_cell_value = graded_complex_ -> value ( top_cell);
      // for (auto x : complex_ -> topstar(v) ) {
      //   if (graded_complex_ -> value (x) > top_cell_value) continue; //MAYBE
      // }
      //}
      if ( maxc != 0 || complex_-> mincoords(v) != 0) continue; //Experimental testing for braids!


      for (Integer type = 0; type < (1L << H); ++ type ) {

        //bound on dimnension? TODO: if maxc = 0, then consider shape as type
        Integer shape = type;
        //if ( maxc == 0) { shape = complex_ -> ST()[type]; }

        //Integer dim = popcount_ ( shape );
        //if ( maxc == 0 and dim > D  ) break; // config space optimization
        //if (dim > D) continue; // config space optimization
        
        //Integer cell = v + complex_ -> type_size() * complex_ -> TS() [ embed(shape, maxc) ];
        //if (graded_complex_ -> value (cell) != 0) continue; // Config space optimization

        if ( ! visited [shape] ) {
          //TODO make local_mate_ an iterative algorithm
          //Provide proof no issue using local_mate_ here then mate_ for determining boundary
          if ( mate_(shape, v, maxc, H, visited) == shape ) {
            //std::cout << " CRIT: " << cell << " shape: " << shape << "\n";
            //std::cout << " CRIT grading: " << graded_complex_ -> value(cell) << "\n";
            Integer dim = popcount_ ( shape ); //Comment out for config
            Integer cell = v + complex_ -> type_size() * complex_ -> TS() [ embed(shape, maxc) ];
            crit_cells [ dim ] . push ( cell );
          }
        }
      }

      // //Doesn't skip fringe
      // Integer H = 1L << D;
      // std::vector<bool> visited ( H );

      // for (Integer type = 0; type < H; ++ type ) { //for each type (hypercube)

      //   Integer shape = complex_ -> ST() [ type ];
      //   Integer c = v + complex_ -> type_size() * type; //cell index

      //   //auto dim = popcount_ ( shape );
      //   //if ( dim > D ) break; //Configuration space optimization

      //   //if (graded_complex_ -> value(c) != 0) continue; //Config space optimization
      //   //std::cout << "c: " << c << " mate: " << mate(c) << "\n";
      //   if ( ! visited [ type ] ) {
      //     if (! complex_ -> rightfringe(c) ) {
      //       //std::cout << " c: " << c << " t: " << type << " mate: " << mate(c) << " visited: " << visited[type] << "\n";
      //       if ( mate_(c,D,visited) == c) { // H. vs D?
      //         auto dim = popcount_ ( shape );
      //         crit_cells[dim] . push ( c );
      //         //std::cout << "c: " << c << " dim: " << dim << "\n";
      //       }
      //     }
      //   }
      // }
    }
    //Unroll critical cells
    Integer idx = 0;
    for (Integer i = 0; i < crit_cells . size(); ++i ) {
      begin_[i] = idx;
      while ( ! crit_cells[i] . empty() ) {
        Integer c = crit_cells[i].front();
        reindex_ . push_back ({c,idx});
        crit_cells[i] . pop();
        ++ idx;
      }
    }
    begin_[D+1] = idx;

    // for ( Integer d = 0; d <= D; ++ d) {
    //   begin_[d] = idx;
    //   for ( auto v : (*complex_)(d) ) { // TODO: skip fringe cells
    //     if ( ! complex_ -> rightfringe(v)  ) { //&& graded_complex_ -> value(v) == 0 
    //       //std::cout << "v: " << v << " mate: " << mate(v) << "\n";
    //       if ( mate(v) == v ) { 
    //         reindex_.push_back({v,idx});
    //         ++idx;
    //       }
    //     }
    //   }
    // }
    // begin_[D+1] = idx;
  }

  /// critical_cells
  std::pair<BeginType const&,ReindexType const&>
  critical_cells ( void ) const {
    return {begin_,reindex_};
  }

  // /// mate
  // Integer
  // mate ( Integer x ) const { 
  //   return mate_(x, complex_ -> dimension());
  // }
  Integer
  mate ( Integer x ) const {
    Integer position = x % complex_ -> type_size(); 
    Integer maxc = complex_ -> maxcoords(position);
    Integer shape = retract(complex_ -> cell_shape(x), maxc);
    Integer H = complex_ -> dimension() - popcount_ ( maxc );
    return position + complex_ -> type_size() * complex_ -> TS() [ embed(mate_(shape,position,maxc,H), maxc) ];
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

  Integer
  popcount_ ( Integer x ) const {
    // http://lemire.me/blog/2016/05/23/the-surprising-cleverness-of-modern-compilers/
    Integer pcnt = 0; 
    while(x != 0) { x &= x - 1; ++pcnt; } 
    return pcnt;
  }

  // Integer
  // embed ( Integer shape, Integer maxc) const {
  //   if (maxc == 0) return shape;

  //   Integer count = 0;
  //   Integer embedding = 0;
  //   Integer embed_coords = ~maxc;
  //   for ( Integer d = 0, bit = 1; d < complex_ -> dimension(); ++ d, bit <<= 1L  ) {
  //     //std::cout << " bit: " << bit << "\n";
  //     if ( embed_coords & bit ) {
  //       if ( shape & (1 << count) ) {
  //         embedding += (1 << d);
  //       }
  //       count += 1;
  //     }
  //   }
  //   return embedding;
  // }
  //TODO example of embed, retract
  //Decode succinct representation (code)
  Integer
  embed ( Integer code, Integer shape) const {
    if (shape == 0) return code;

    Integer count = 0;
    Integer stub = 0;
    Integer temp = shape;
    while (code != 0) {
      if (!(shape & 1)) {
        if (code & 1) stub += (1<<count);
        code >>= 1;
      }
      shape >>= 1;
      ++count; 
    }
    return stub;
  }
  //Retract function, e.g., retract(embed(code,shape),shape) = code
  //Create succinct representation
  Integer
  retract (Integer code, Integer shape ) const {
    if (shape == 0) return code;

    Integer count = 0;
    Integer stub = 0;
    while (code != 0) {
      if (!(shape & 1)) {
        if (code & 1) stub += (1<<count);
        ++count;
      }
      code >>= 1;
      shape >>= 1;
    }
    return stub;
  }

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
    
    //This does not seem to be needed 
    if ( complex_ -> rightfringe(cell) ) { 
      //std::cout << "Fringe in mate_\n";
      return cell; // MAYBE
    }

    //Integer mincoords = complex_ -> mincoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
    //Integer maxcoords = complex_ -> maxcoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
    Integer shape = complex_ -> cell_shape(cell);
    Integer position = cell % complex_ -> type_size();
    if ( position == complex_ -> type_size() - 1 ) return cell; // Break cycles
    for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= 1L  ) {
      // If on right fringe for last dimension, prevent mating with left boundary
      if ( (d == D-1) && (position + complex_ -> PV()[d] >= complex_ -> type_size()) ) {
        //std::cout << "Position in mate_\n";
        break; //Does this ever run?
      }
      //if ( fringe && (mincoords & bit) ) continue; // Todo -- is this the best
      //if ( bit & maxcoords ) continue; // Don't connect fringe to acyclic part
      Integer type_offset = complex_ -> type_size() * complex_ -> TS() [ shape ^ bit ];
      Integer proposed_mate = position + type_offset;
      //std::cout << "position : " << position << " type_offset: " << type_offset << "\n";
      if ( graded_complex_ -> value(proposed_mate) == graded_complex_ -> value(cell) && proposed_mate == mate_(proposed_mate, d) ) { 
        return proposed_mate;
      }
    }
    return cell;
  }

  // Integer mate_ ( Integer cell, Integer D, std::vector<bool> & visited ) const {
  //   //bool fringe = complex_ -> rightfringe(cell);
    
  //   //This does not seem to be needed if mate is not run on fringe
  //   // if ( complex_ -> rightfringe(cell) ) return cell; // MAYBE
  //   Integer type = complex_ -> cell_type ( cell );
  //   Integer shape = complex_ -> ST() [ type ];

  //   visited[type] = true; //mark visited
  //   //Integer mincoords = complex_ -> mincoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
  //   //Integer maxcoords = complex_ -> maxcoords(cell); // TODO: optimize to compute this as it loops through d rather than demanding all
  //   //Integer shape = complex_ -> cell_shape(cell);
  //   Integer position = cell % complex_ -> type_size();
  //   if ( position == complex_ -> type_size() - 1 ) return cell; // Break cycles
  //   for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= 1L  ) {
  //     // If on right fringe for last dimension, prevent mating with left boundary
  //     if ( (d == D-1) && (position + complex_ -> PV()[d] >= complex_ -> type_size()) )  break; //Does this ever run?
  //     //if ( fringe && (mincoords & bit) ) continue; // Todo -- is this the best
  //     //if ( bit & maxcoords ) continue; // Don't connect fringe to acyclic part
  //     Integer type_offset = complex_ -> type_size() * complex_ -> TS() [ shape ^ bit ];
  //     Integer proposed_mate = position + type_offset;
  //     //std::cout << "position : " << position << " type_offset: " << type_offset << "\n";
  //     if ( graded_complex_ -> value(proposed_mate) == graded_complex_ -> value(cell) && proposed_mate == mate_(proposed_mate, d, visited) ) { 
  //       return proposed_mate;
  //     }
  //   }
  //   return cell;
  // }
  Integer mate_ ( Integer shape, Integer position, Integer maxc, Integer D) const {

    if ( position == complex_ -> type_size() - 1 ) return shape; // Break cycles
    for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= 1L  ) {
      // If on right fringe for last dimension, prevent mating with left boundary
      //if ( (d == D-1) && (position + complex_ -> PV()[d] >= complex_ -> type_size()) )  break; //Does this ever run?
      Integer shape_offset = shape ^ bit;

      Integer cell = position + complex_ -> type_size() * complex_ -> TS() [ embed(shape, maxc) ];
      Integer proposed_mate = position + complex_ -> type_size() * complex_ -> TS() [ embed(shape_offset, maxc) ];
      //std::cout << " mate_cell: " << cell << " mate_proposed_mate: " << proposed_mate << "\n";
      //std::cout << " v(mate): " << graded_complex_ -> value(proposed_mate) << " v(cell): " << graded_complex_ -> value(cell) << "\n";
      if ( graded_complex_ -> value(proposed_mate) == graded_complex_ -> value(cell) && shape_offset == mate_(shape_offset, position, maxc, d) ) { 
        return shape_offset;
      }
    }
    return shape;
  }
  //DRY SMELL
  Integer mate_ ( Integer shape, Integer position, Integer maxc, Integer D, std::vector<bool> & visited ) const {
    
    visited[shape] = true; //mark visited
    if ( position == complex_ -> type_size() - 1 ) return shape; // Break cycles
    for ( Integer d = 0, bit = 1; d < D; ++ d, bit <<= 1L  ) {
      // If on right fringe for last dimension, prevent mating with left boundary
      //if ( (d == D-1) && (position + complex_ -> PV()[d] >= complex_ -> type_size()) )  break; //Does this ever run?
      Integer shape_offset = shape ^ bit;

      Integer cell = position + complex_ -> type_size() * complex_ -> TS() [ embed(shape, maxc) ];
      Integer proposed_mate = position + complex_ -> type_size() * complex_ -> TS() [ embed(shape_offset, maxc) ];
      //std::cout << " cell: " << cell << " proposed_mate: " << proposed_mate << "\n";
      //std::cout << " v(mate): " << graded_complex_ -> value(proposed_mate) << " v(cell): " << graded_complex_ -> value(cell) << "\n";
      if ( graded_complex_ -> value(proposed_mate) == graded_complex_ -> value(cell) && shape_offset == mate_(shape_offset, position, maxc, d, visited) ) { 
        return shape_offset;
      }
    }
    return shape;
  }
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void
CubicalNNMorseMatchingBinding(py::module &m) {
  py::class_<CubicalNNMorseMatching, std::shared_ptr<CubicalNNMorseMatching>>(m, "CubicalNNMorseMatching")
    .def(py::init<std::shared_ptr<CubicalComplex>>())
    .def(py::init<std::shared_ptr<GradedComplex>>())    
    .def("mate", &CubicalNNMorseMatching::mate)
    .def("priority", &CubicalNNMorseMatching::priority);
}
