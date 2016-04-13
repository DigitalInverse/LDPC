#ifndef __ldpc_802_11_utilities_h
#define __ldpc_802_11_utilities_h

#include <iostream>

// Piecewise linear approximation of tanh which
// mitigates error floor caused by std::tanh.
//
const double tanh_lin(const double x);

// Quantization of tanh
//
const double tanh_quant(const double x);

// Piecewise linear approximation of atanh which
// mitigates error floor caused by std::atanh.
//
const double atanh_lin(const double x);

// Quantization of atanh
//
const double atanh_quant(const double x);

// Somewhat nicer formatting than itpp::cout for mat
//
template<class T>
void pretty_print_matrix(const T & M)
{
  for (int row=0; row<M.rows(); ++row) {
    for (int col=0; col<M.cols(); ++col) {
      std::cout << M(row,col) << "\t";
    }
    std::cout << std::endl;
  }
}

#endif //__ldpc_802_11_utilities_h
