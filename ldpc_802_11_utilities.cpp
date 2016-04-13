#include "ldpc_802_11_utilities.h"

// Piecewise linear approximation of tanh which
// avoids error floor caused by std::tanh. The
// values are from
//
//    S. Papaharalabos, P. Sweeney, B.G. Evans,
//    P.T. Mathiopoulos, G. Albertazzi, A. Vanelli-Coralli,
//    and G.E. Corazza.
//    "Modified Sum-Product Algorithms for Decoding
//    Low-Density Parity-Check Codes,"
//    IET Communications, vol. 1, no. 3, pp. 294--300,
//    June 2007.
//
// The implementation is computationally faster than
// std::tanh (exploiting the symmetry of tanh does
// not speedup the following naive implementation).
//
const double tanh_lin(const double x)
{
  if (x <= -7.0) {
    return -0.999998;
  }
  else if (x > 7.0) {
    return 0.999998;
  }
  else if (x <= -3.0) {
    return 0.0012*x - 0.9914;
  }
  else if (x <= -1.6) {
    return 0.0524*x - 0.8378;
  }
  else if (x <= -0.8) {
    return 0.322*x - 0.4064;
  }
  else if (x <= 0.8) {
    return 0.83*x;
  }
  else if (x <= 1.6) {
    return 0.322*x + 0.4064;
  }
  else if (x <= 3.0) {
    return 0.0524*x + 0.8378;
  }
  else { // if (x <= 7.0) {
    return 0.0012*x + 0.9914;
  }
}

// Quantization of tanh, as described in
//
//    S. Papaharalabos, P. Sweeney, B.G. Evans,
//    P.T. Mathiopoulos, G. Albertazzi, A. Vanelli-Coralli,
//    and G.E. Corazza.
//    "Modified Sum-Product Algorithms for Decoding
//    Low-Density Parity-Check Codes,"
//    IET Communications, vol. 1, no. 3, pp. 294--300,
//    June 2007.
//
// Not as good as linear approximation, and
// somewhat slower (possibly because one of
// more comparison).
//
const double tanh_quant(const double x)
{
  if (x <= -7.0) {
    return -0.999998;
  }
  else if (x > 7.0) {
    return 0.999998;
  }
  else if (x <= -3.0) {
    return -0.99991;
  }
  else if (x <= -1.6) {
    return -0.9801;
  }
  else if (x <= -0.8) {
    return -0.8337;
  }
  else if (x <= 0) {
    return -0.3799;
  }
  else if (x <= 0.8) {
    return 0.3799;
  }
  else if (x <= 1.6) {
    return 0.8337;
  }
  else if (x <= 3.0) {
    return 0.9801;
  }
  else { // if (x <= 7.0) {
    return 0.99991;
  }
}

// Piecewise linear approximation of atanh which
// avoids error floor caused by std::atanh. The
// values are from
//
//    S. Papaharalabos, P. Sweeney, B.G. Evans,
//    P.T. Mathiopoulos, G. Albertazzi, A. Vanelli-Coralli,
//    and G.E. Corazza.
//    "Modified Sum-Product Algorithms for Decoding
//    Low-Density Parity-Check Codes,"
//    IET Communications, vol. 1, no. 3, pp. 294--300,
//    June 2007.
//
// Computationally faster than std::atanh
// (exploiting the symmetry of atanh does not
// speedup the following naive implementation).
//
const double atanh_lin(const double x)
{
  if (x <= -0.999998) {
    return -7.0;
  }
  else if (x > 0.999998) {
    return 7.0;
  }
  else if (x <= -0.9951) {
    return (x + 0.9914)/0.0012;
  }
  else if (x <= -0.9217) {
    return (x + 0.8378)/0.0524;
  }
  else if (x <= -0.6640) {
    return (x + 0.4064)/0.322;
  }
  else if (x <= 0.6640) {
    return x/0.83;
  }
  else if (x <= 0.9217) {
    return (x - 0.4064)/0.322;
  }
  else if (x <= 0.9951) {
    return (x - 0.8378)/0.0524;
  }
  else { //if (x <= 0.999998) {
    return (x - 0.9914)/0.0012;
  }
}

// Quantization of atanh, as described in
//
//    S. Papaharalabos, P. Sweeney, B.G. Evans,
//    P.T. Mathiopoulos, G. Albertazzi, A. Vanelli-Coralli,
//    and G.E. Corazza.
//    "Modified Sum-Product Algorithms for Decoding
//    Low-Density Parity-Check Codes,"
//    IET Communications, vol. 1, no. 3, pp. 294--300,
//    June 2007.
//
// Not as good as linear approximation, and
// somewhat slower (possibly because one of
// more comparison).
//
const double atanh_quant(const double x)
{
  if (x <= -0.999998) {
    return -7.0;
  }
  else if (x > 0.999998) {
    return 7.0;
  }
  else if (x <= -0.9951) {
    return -3.3516;
  }
  else if (x <= -0.9217) {
    return -1.9259;
  }
  else if (x <= -0.6640) {
    return -1.0791;
  }
  else if (x <= 0) {
    return -0.3451;
  }
  else if (x <= 0.6640) {
    return 0.3451;
  }
  else if (x <= 0.9217) {
    return 1.0791;
  }
  else if (x <= 0.9951) {
    return 1.9259;
  }
  else { //if (x <= 0.999998) {
    return 3.3516;
  }
}
