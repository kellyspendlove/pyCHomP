/// ValuationConfig.h
/// Kelly Spendlove
/// 2018-10-1
/// MIT LICENSE

// A specific valuation for calculating the configuration space
// of n unit hard squares in a pxq rectangle

#pragma once

#include "common.h"
//#include <unordered_set>
#include <array>

// std::function<Integer(Integer)>
// construct_config_valuation ( std::shared_ptr<CubicalComplex> c ) {
//   return [=](Integer cell) {
//       if ( c -> rightfringe(cell) ) {
//           return 2;
//       }
//       auto x = c -> barycenter(cell);
//       for ( int i = 0; i < x.size(); i=i+2 ) {
//           int x1 = x[i];
//           int y1 = x[i+1];
//           for ( int j = i+2; j < x.size(); j=j+2) {
//             int x2 = x[j];
//             int y2 = x[j+1];
//             //if both have extent
//             if ( x1 % 2 == 1 and x2 % 2 == 1 ) {
//                 if ( abs(x1-x2) >= 4 ) continue;
//             }
//             else if (x1 % 2 == 1 or x2 % 2 == 1) {
//                 if ( abs(x1-x2) >= 3 ) continue;
//             }
//             else {
//                 if ( abs(x1-x2) >= 2 ) continue;
//             }
//             //if both have extent
//             if ( y1 % 2 == 1 and y2 % 2 == 1 ) {
//                 if ( abs(y1-y2) >= 4 ) continue;
//             }
//             else if (y1 % 2 == 1 or y2 % 2 == 1) {
//                 if ( abs(y1-y2) >= 3 ) continue;
//             }
//             else {
//                 if ( abs(y1-y2) >= 2 ) continue;
//             }
//             return 1;
//           }
//       }
//       // for (int i = 0; i < x.size(); i=i+2) {
//       //   if (x[i+1] > 6) {
//       //     return 1;
//       //   }
//       // }
//       return 0;
//   };
// }
std::function<Integer(Integer)>
construct_config_valuation_order ( std::shared_ptr<CubicalComplex> c ) {
  return [=](Integer cell) {
      // if ( c -> rightfringe(cell) ) { //Can this be removed?
      //   return 2;
      // }
      auto x = c -> barycenter(cell);
      int N = x.size() / 2;
      for ( int i = 0; i < N; ++i ) {
          int x1 = x[i];
          int y1 = x[i+N];
          for ( int j = i+1; j < N; ++j ) {
            int x2 = x[j];
            int y2 = x[j+N];
            //if both have extent
            if ( x1 % 2 == 1 and x2 % 2 == 1 ) {
                if ( abs(x1-x2) >= 4 ) continue;
            }
            else if (x1 % 2 == 1 or x2 % 2 == 1) {
                if ( abs(x1-x2) >= 3 ) continue;
            }
            else {
                if ( abs(x1-x2) >= 2 ) continue;
            }
            //if both have extent
            if ( y1 % 2 == 1 and y2 % 2 == 1 ) {
                if ( abs(y1-y2) >= 4 ) continue;
            }
            else if (y1 % 2 == 1 or y2 % 2 == 1) {
                if ( abs(y1-y2) >= 3 ) continue;
            }
            else {
                if ( abs(y1-y2) >= 2 ) continue;
            }
            return 1;
          }
      }
      return 0;
  };
}
std::function<Integer(Integer)>
construct_config_valuation ( std::shared_ptr<CubicalComplex> c ) {
  return [=](Integer cell) {
    
    auto result = c -> coordinates(cell);
    auto shape = c -> cell_shape ( cell );
    Integer N = result . size() / 2;
    Integer p = c -> boxes()[0];
    Integer q = c -> boxes()[N];
    auto shape_y = shape >> N;
    //std::vector<bool> board ( p*q ); //positions of squares
    std::array<bool,128> board {}; //positions of squares (up to 11x11 board)
    for (Integer i = 0, bit_x = shape, bit_y = shape_y; i < N; ++ i, bit_x >>= 1, bit_y >>= 1 ) {
      Integer x = result[i];
      Integer y = result[i+N];
      //Extents
      // Integer e_x = (shape >> i) & 1;
      // Integer e_y = (shape >> (i+N)) & 1;
      Integer e_x = bit_x & 1;
      Integer e_y = bit_y & 1;
      if ( board [x+(p+1)*y] ) return 1;
      if ( e_x && board [(x+e_x)+(p+1)*y] ) return 1;
      if ( e_y && board [x + (p+1)*(y+e_y)] ) return 1;
      if ( e_x && e_y && board [(x+e_x) + (p+1)*(y+e_y)] ) return 1;

      //std::cout << " no collision \n";
      board [x+(p+1)*y] = true;
      if ( e_x ) board [(x+e_x)+(p+1)*y ] = true;
      if ( e_y ) board [x + (p+1)*(y+e_y)] = true;
      if ( e_x && e_y) board [ (x+e_x) + (p+1)*(y+e_y)] = true;

    }
    return 0;
  };
}
std::function<Integer(Integer)>
grading_config_cubes ( std::shared_ptr<CubicalComplex> c ) {
  return [=](Integer cell) {
      if ( c -> rightfringe(cell) ) {
          return 2;
      }
      auto x = c -> barycenter(cell);
      int N = x.size() / 3;
      for ( int i = 0; i < N; ++i ) {
          int x1 = x[i];
          int y1 = x[i+N];
          int z1 = x[i+2*N];
          for ( int j = i+1; j < N; ++j ) {
            int x2 = x[j];
            int y2 = x[j+N];
            int z2 = x[j+2*N];
            //x checks
            if ( x1 % 2 == 1 and x2 % 2 == 1 ) {
                if ( abs(x1-x2) >= 4 ) continue;
            }
            else if (x1 % 2 == 1 or x2 % 2 == 1) {
                if ( abs(x1-x2) >= 3 ) continue;
            }
            else {
                if ( abs(x1-x2) >= 2 ) continue;
            }
            //y checks
            if ( y1 % 2 == 1 and y2 % 2 == 1 ) {
                if ( abs(y1-y2) >= 4 ) continue;
            }
            else if (y1 % 2 == 1 or y2 % 2 == 1) {
                if ( abs(y1-y2) >= 3 ) continue;
            }
            else {
                if ( abs(y1-y2) >= 2 ) continue;
            }
            //z checks
            if ( z1 % 2 == 1 and z2 % 2 == 1 ) {
                if ( abs(z1-z2) >= 4 ) continue;
            }
            else if (z1 % 2 == 1 or z2 % 2 == 1) {
                if ( abs(z1-z2) >= 3 ) continue;
            }
            else {
                if ( abs(z1-z2) >= 2 ) continue;
            }
            return 1;
          }
      }
      return 0;
  };
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

inline void
ValuationConfigBinding(py::module &m) {
  m.def("construct_config_valuation", &construct_config_valuation);
  m.def("construct_config_valuation_order", &construct_config_valuation_order);
  m.def("grading_config_cubes", &grading_config_cubes);
}

