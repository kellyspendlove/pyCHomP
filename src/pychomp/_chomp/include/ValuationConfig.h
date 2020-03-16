/// ValuationConfig.h
/// Kelly Spendlove
/// 2018-10-1
/// MIT LICENSE

// A specific valuation for calculating the configuration space
// of n unit hard squares in a pxq rectangle

#pragma once

#include "common.h"

std::function<Integer(Integer)>
construct_config_valuation ( std::shared_ptr<CubicalComplex> c ) {
  return [=](Integer cell) {
      if ( c -> rightfringe(cell) ) {
          return 2;
      }
      auto x = c -> barycenter(cell);
      for ( int i = 0; i < x.size(); i=i+2 ) {
          int x1 = x[i];
          int y1 = x[i+1];
          for ( int j = i+2; j < x.size(); j=j+2) {
            int x2 = x[j];
            int y2 = x[j+1];
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
      // for (int i = 0; i < x.size(); i=i+2) {
      //   if (x[i+1] > 6) {
      //     return 1;
      //   }
      // }
      return 0;
  };
}
std::function<Integer(Integer)>
construct_config_valuation_order ( std::shared_ptr<CubicalComplex> c ) {
  return [=](Integer cell) {
      if ( c -> rightfringe(cell) ) {
          return 2;
      }
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

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

inline void
ValuationConfigBinding(py::module &m) {
  m.def("construct_config_valuation", &construct_config_valuation);
  m.def("construct_config_valuation_order", &construct_config_valuation_order);
}

