#include "ldpc_802_11.h"
#include "ldpc_802_11_utilities.h"

// constructor which creates a parity check matrix
// and the Tanner graph representation of the code
//
LDPC_802_11_codec::LDPC_802_11_codec(const int block_length,
                                     const code_rate_type rate,
                                     const decode_method_type decode_method,
                                     const bool debug) :
          _block_length(block_length),
          _rate(rate),
          _decode_method(decode_method)
{
  _generate_matrix_prototype(block_length,
                             rate,
                             debug);
  _compute_row_x(debug);
  _create_parity_check_matrix(debug);
  _create_tanner_graph(debug);
}


// The encoding algorithm that is implemented below is
// described in
//
//    Z. Cai, J. Hao, P.H. Tan, S. Sun, and P.S. Chin,
//    "Efficient encoding of IEEE 802.11n LDPC codes,"
//    Electronics Letters, vol. 42, no. 25, pp. 1471--1472,
//    Dec. 2006.
//
// This scheme is approximately 23 times faster than naive encoding
// (multiplying the vector of information bits with the generator
// matrix) for a block size of 1944 bits (the biggest block length)
// and approximately 5 times faster for a block length of 648 bits
// (the smallest block length).
//
void LDPC_802_11_codec::encode(itpp::bvec & output_bits,
                               const itpp::bvec & input_bits) const
{
  output_bits = itpp::zeros_b(_N);
  output_bits.set_subvector(0, input_bits);
  int position = _k_b * _Z;

  // compute p_0
  itpp::bvec p_0 = itpp::zeros_b(_Z);
  for (int i=0; i<_m_b; ++i) {
    p_0 += _compute_lambda_i(i, input_bits);
  }
  output_bits.set_subvector(position, p_0);
  position += _Z;

  // compute p_1
  itpp::bvec shifted_p_0 = p_0;
  shifted_p_0.shift_left(shifted_p_0.left(1)); // \Pi_1 p_0 in the paper
  itpp::bvec p_1 = _compute_lambda_i(0, input_bits) + shifted_p_0;
  output_bits.set_subvector(position, p_1);
  position += _Z;

  // upward recursion p_{i+1} = p_i + \lambda_i
  itpp::bvec p_i, p_i_1; // p_i and p_{i+1}
  p_i = p_1;
  for (int i=1; i<_row_x; ++i) {
    p_i_1 = p_i + _compute_lambda_i(i, input_bits);
    output_bits.set_subvector(position, p_i_1);
    position += _Z;
    p_i = p_i_1;
  }

  // compute p_{m_b-1}
  position = (_n_b -1) * _Z;
  itpp::bvec p_last = _compute_lambda_i(_m_b-1, input_bits) + shifted_p_0; // p_{m_b-1}
  output_bits.set_subvector(position, p_last);
  position -= _Z;

  // downward recursion p_i = p_{i+1} + \lambda_i
  p_i_1 = p_last;
  for (int i=_m_b-2; i>_row_x; --i) {
    p_i = p_i_1 + _compute_lambda_i(i, input_bits);
    output_bits.set_subvector(position, p_i);
    position -= _Z;
    p_i_1 = p_i;
  }

  // row x
  p_i = _compute_lambda_i(_row_x, input_bits) + p_0 + p_i_1;
  output_bits.set_subvector(position, p_i);
}

// set reliability factor that is used for soft scaling of the input
//
void LDPC_802_11_codec::set_channel_reliability_value(const double noise_variance)
{
  _L_c = 2.0 / noise_variance;
}

// helper function used in the encoding method
//
itpp::bvec LDPC_802_11_codec::_compute_lambda_i(const int i,
                                                const itpp::bvec m) const
{
  itpp::bvec lambda_i = itpp::zeros_b(_Z);
  for (int j=0; j<_k_b; ++j) {
    itpp::bvec m_j = m.get(j*_Z, (j+1)*_Z-1);
    int shift = _matrix_prototype(i, j);
    if (shift != -1) {
      if (shift != 0) {
        m_j.shift_left(m_j.left(shift)); // circular shift of m_j
      }
      lambda_i += m_j;
    }
  }
  return lambda_i;
}


// Compute parity check matrix prototype specified
// in 802.11-2012 Annex F, pp. 2304--2306,
// http://standards.ieee.org/getieee802/download/802.11-2012.pdf
//
void LDPC_802_11_codec::_generate_matrix_prototype(const int block_length = 648,
                                                   const code_rate_type rate = half,
                                                   const bool print = false)
{
  if (block_length == 648) { // -----------------------------------------

    _Z = 27;
    _m_b = 12;
    _n_b = 24;
    _k_b = _n_b - _m_b;

    _matrix_prototype = -1 * itpp::ones(_m_b, _n_b);

    switch (rate) {
    case half: {
      _matrix_prototype.set_col(0, "0 22 6 2 23 24 25 13 7 11 25 3");
      _matrix_prototype.set(1,1, 0);
      _matrix_prototype.set(7,1, 24);
      _matrix_prototype.set(8,1, 20);
      _matrix_prototype.set(2,2, 0);
      _matrix_prototype.set(5,2, 23);
      _matrix_prototype.set(10,2, 8);
      _matrix_prototype.set(3,3, 0);
      _matrix_prototype.set(5,3, 1);
      _matrix_prototype.set(8,3, 16);
      _matrix_prototype.set_col(4, "0 17 10 20 3 17 8 0 22 19 23 16");
      _matrix_prototype.set(0,5, 0);
      _matrix_prototype.set(8,5, 10);
      _matrix_prototype.set(10,5, 18);
      _matrix_prototype.set(1,6, 0);
      _matrix_prototype.set(5,6, 3);
      _matrix_prototype.set(7,6, 8);
      _matrix_prototype.set(1,7, 0);
      _matrix_prototype.set(10,7, 14);
      _matrix_prototype.set(11,7, 2);
      _matrix_prototype.set_col(8, "0 12 24 25 0 10 7 6 23 13 9 25");
      _matrix_prototype.set(3,9, 0);
      _matrix_prototype.set(6,9, 18);
      _matrix_prototype.set(11,9, 5);
      _matrix_prototype.set(2,10, 0);
      _matrix_prototype.set(4,10, 9);
      _matrix_prototype.set(9,10, 3);
      _matrix_prototype.set(0,11, 0);
      _matrix_prototype.set(4,11, 11);
      _matrix_prototype.set(9,11, 17);
      _matrix_prototype.set(0,12, 1);
      _matrix_prototype.set(6,12, 0);
      _matrix_prototype.set(11,12, 1);
      break;
    }
    case two_thirds:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    case three_quarters:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    case five_sixths:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    default:
      it_error("No such LDPC code rate in IEEE 802.11!");
      break;
    }
  }
  else if (block_length == 1296) { // -----------------------------------

    _Z = 54;
    _m_b = 12;
    _n_b = 24;
    _k_b = _n_b - _m_b;

    _matrix_prototype = -1 * itpp::ones(_m_b, _n_b);

    switch (rate) {
    case half:
      _matrix_prototype.set_col(0, "40 50 39 33 45 51 47 5 33 1 -1 49");
      _matrix_prototype.set(1,1, 1);
      _matrix_prototype.set(2,1, 50);
      _matrix_prototype.set(6,1, 11);
      _matrix_prototype.set(10,1, 18);
      _matrix_prototype.set(7,2, 25);
      _matrix_prototype.set(9,2, 27);
      _matrix_prototype.set(11,2, 17);
      _matrix_prototype.set(3,3, 38);
      _matrix_prototype.set(5,3, 48);
      _matrix_prototype.set(8,3, 34);
      _matrix_prototype.set_col(4, "22 48 4 37 0 35 -1 6 24 1 23 30");
      _matrix_prototype.set(1,5, 35);
      _matrix_prototype.set(4,5, 22);
      _matrix_prototype.set(6,5, 17);
      _matrix_prototype.set(0,6, 49);
      _matrix_prototype.set(2,6, 2);
      _matrix_prototype.set(7,6, 45);
      _matrix_prototype.set(0,7, 23);
      _matrix_prototype.set(3,7, 4);
      _matrix_prototype.set(10,7, 8);
      _matrix_prototype.set_col(8, "43 13 -1 1 20 44 51 13 23 38 0 34");
      _matrix_prototype.set(4,9, 42);
      _matrix_prototype.set(7,9, 40);
      _matrix_prototype.set(10,9, 35);
      _matrix_prototype.set(1,10, 30);
      _matrix_prototype.set(5,10, 18);
      _matrix_prototype.set(9,10, 44);
      _matrix_prototype.set(2,11, 49);
      _matrix_prototype.set(8,11, 46);
      _matrix_prototype.set(11,11, 19);
      _matrix_prototype.set(0,12, 1);
      _matrix_prototype.set(6,12, 0);
      _matrix_prototype.set(11,12, 1);
      break;
    case two_thirds:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    case three_quarters:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    case five_sixths:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    default:
      it_error("No such LDPC code rate in IEEE 802.11!");
      break;
    }
  }
  else if (block_length == 1944) { // -----------------------------------

    _Z = 81;
    _m_b = 12;
    _n_b = 24;
    _k_b = _n_b - _m_b;

    _matrix_prototype = -1 * itpp::ones(_m_b, _n_b);

    switch (rate) {
    case half:
      _matrix_prototype.set_col(0, "57 3 30 62 40 0 69 65 64 -1 2 24");
      _matrix_prototype.set(3,1, 53);
      _matrix_prototype.set(6,1, 79);
      _matrix_prototype.set(9,1, 45);
      _matrix_prototype.set(10,1, 56);
      _matrix_prototype.set(1,2, 28);
      _matrix_prototype.set(6,2, 79);
      _matrix_prototype.set(11,2, 61);
      _matrix_prototype.set(4,3, 20);
      _matrix_prototype.set(9,3, 70);
      _matrix_prototype.set(10,3, 57);
      _matrix_prototype.set_col(4, "50 0 24 53 66 8 -1 38 14 0 35 60");
      _matrix_prototype.set(2,5, 37);
      _matrix_prototype.set(7,5, 57);
      _matrix_prototype.set(8,5, 52);
      _matrix_prototype.set(0,6, 11);
      _matrix_prototype.set(5,6, 42);
      _matrix_prototype.set(6,6, 56);
      _matrix_prototype.set(3,7, 3);
      _matrix_prototype.set(4,7, 22);
      _matrix_prototype.set(11,7, 27);
      _matrix_prototype.set_col(8, "50 55 56 35 28 50 52 72 30 77 -1 51");
      _matrix_prototype.set(1,9, 7);
      _matrix_prototype.set(2,9, 14);
      _matrix_prototype.set(9,9, 9);
      _matrix_prototype.set(0,10, 79);
      _matrix_prototype.set(7,10, 27);
      _matrix_prototype.set(10,10, 12);
      _matrix_prototype.set(5,11, 8);
      _matrix_prototype.set(8,11, 32);
      _matrix_prototype.set(11,11, 16);
      _matrix_prototype.set(0,12, 1);
      _matrix_prototype.set(6,12, 0);
      _matrix_prototype.set(11,12, 1);
      break;
    case two_thirds:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    case three_quarters:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    case five_sixths:
      // TODO Implement this rate (see IEEE standard)
      it_error("Implement this!");
      break;
    default:
      it_error("No such LDPC code rate in IEEE 802.11!");
      break;
    }
  }
  else {
    it_error("No such LDPC block length IEEE 802.11!");
  }

  // compute rightmost part of parity check matrix prototype
  int row_ctr = 0;
  for (int col_ctr=_k_b+1; col_ctr<_n_b; ++col_ctr) {
    _matrix_prototype.set(row_ctr,col_ctr, 0);
    row_ctr++;
    _matrix_prototype.set(row_ctr,col_ctr, 0);
  }

  _M = _m_b * _Z;
  _N = _n_b * _Z;

  if (print) {
    std::cout << "parity check matrix prototype:\n";
    pretty_print_matrix(_matrix_prototype);
  }
}


// helper function to find _row_x, which is used in the
// encoding method (see paper about encoding listed above)
//
void LDPC_802_11_codec::_compute_row_x(const bool print_row_x = false)
{
  for (int row=1; row<_m_b; ++row) {
    if (_matrix_prototype(row,_k_b)==0) {
      _row_x = row;
    }
  }
  if (print_row_x) {
    std::cout << "row x: " << _row_x << std::endl;
  }
}


// build up the parity check matrix, which, in turn, is used
// to create the Tanner graph representation, as well as to
// verify if a decoded bit vector is a valid codeword or not
//
void LDPC_802_11_codec::_create_parity_check_matrix(const bool print = false)
{
  _parity_check_matrix = itpp::zeros_b(_M, _N);

  // compute leftmost part of the parity check matrix
  itpp::bmat eye_Z = itpp::eye_b(_Z);
  itpp::bmat eye_Z_shifted(_Z, _Z);

  for (int row=0; row<_m_b; ++row) {
    for (int col=0; col<_n_b; ++col) {
      int shift = _matrix_prototype(row, col);
      if (shift != -1) {
        if (shift == 0) {
          _parity_check_matrix.set_submatrix(row*_Z, col*_Z, eye_Z);
        }
        else {
          eye_Z_shifted.set_submatrix(0, 0, eye_Z.get(0, _Z-1, _Z - shift, _Z-1));    // leftmost part
          eye_Z_shifted.set_submatrix(0, shift, eye_Z.get(0, _Z-1, 0, _Z-1 - shift)); // rightmost part
          _parity_check_matrix.set_submatrix(row*_Z, col*_Z, eye_Z_shifted);
        }
      }
    } // col
  } // row

  if (print) {
    std::cout << "parity check matrix:\n";
    pretty_print_matrix(_parity_check_matrix);
  }
}


// create the Tanner graph representation of the code
// (a bipartite graph with check nodes and symbol nodes)
//
void LDPC_802_11_codec::_create_tanner_graph(const bool print = false)
{
  // build a set of check nodes, where each check node
  // is connected to a set of symbol/variable nodes
  _check_node_array.set_size(_parity_check_matrix.rows());

  for (int row=0; row<_parity_check_matrix.rows(); ++row) {
    itpp::ivec tmp = ""; // empty vector
    for (int col=0; col<_parity_check_matrix.cols(); ++col) {
      if (_parity_check_matrix(row, col) == 1) {
        tmp = itpp::concat(tmp, col);
      }
    }
    _check_node_array(row) = tmp;
  }

  // build a set of symbol nodes, where each symbol node
  // is connected to a set of check nodes
  _symbol_node_array.set_size(_parity_check_matrix.cols());

  for (int col=0; col<_parity_check_matrix.cols(); ++col) {
    itpp::ivec tmp = ""; // empty vector
    for (int row=0; row<_parity_check_matrix.rows(); ++row) {
      if (_parity_check_matrix(row, col) == 1) {
        tmp = itpp::concat(tmp, row);
      }
    }
    _symbol_node_array(col) = tmp;
  }

  if (print) {
    std::cout << _check_node_array << std::endl;
    std::cout << _symbol_node_array << std::endl;
  }
}


int LDPC_802_11_codec::get_block_length() const
{
  return _block_length;
}


void LDPC_802_11_codec::get_code_rate(int & k, int & n) const
{
  switch (_rate) {
  case half:
    k = 1;
    n = 2;
    break;
  case two_thirds:
    k = 2;
    n = 3;
    break;
  case three_quarters:
    k = 3;
    n = 4;
    break;
  case five_sixths:
    k = 5;
    n = 6;
    break;
  default:
    it_error("No such LDPC code rate in IEEE 802.11!");
    break;
  }
}


// check if a given bit vector is a valid codeword or not
//
bool LDPC_802_11_codec::verify_codeword(const itpp::bvec & codeword) const
{
  bool is_codeword = false;

  itpp::bmat cw_bmat = itpp::zeros_b(1, _N);
  cw_bmat.set_row(0, codeword);

  if ((cw_bmat * _parity_check_matrix.transpose()) == itpp::zeros_b(1, _M)) {
    is_codeword = true;
  }

  return is_codeword;
}


// Message passing decoding of an LDPC code based on the
// flooding scheme. We follow the notation in
//
//    S. Papaharalabos, P. Sweeney, B.G. Evans,
//    P.T. Mathiopoulos, G. Albertazzi, A. Vanelli-Coralli,
//    and G.E. Corazza.
//    "Modified Sum-Product Algorithms for Decoding
//    Low-Density Parity-Check Codes,"
//    IET Communications, vol. 1, no. 3, pp. 294--300,
//    June 2007.
//
void LDPC_802_11_codec::decode_flooding(itpp::bvec & decoded_bits,
                                        const itpp::vec & received_symbols,
                                        const int max_iterations,
                                        const int check_iterations)
{
  int stop_criterion = 0;                  // start to check
  if (max_iterations > check_iterations) { // if valid codeword
    stop_criterion = check_iterations;     // after stop_criterion
  }                                        // iterations

  // 1. Initialization step:
  //
  itpp::mat lambda_n_m = itpp::zeros(_M, _N);
  itpp::mat Lambda_m_n = itpp::zeros(_M, _N);

  itpp::vec L_u_n = _L_c * received_symbols;
  for (int n=0; n<_N; ++n) {
    itpp::vec tmp = L_u_n(n) * itpp::ones(_M);
    lambda_n_m.set_col(n, tmp);
  }

  // 2. Iteration step:
  //
  int iter = 1;
  bool corrected = false;
  while (!(corrected) && (iter<=max_iterations)) {

    // Part A: Message passing from check nodes to symbol nodes
    for (int m=0; m<_M; ++m) {
      itpp::ivec N_m = _check_node_array(m); // set of symbol nodes connected to check node m
      for (int n_idx=0; n_idx<N_m.size(); ++n_idx) {
        if (_decode_method == product_version) {
          Lambda_m_n(m, N_m(n_idx)) = _product_version(N_m, lambda_n_m, m, n_idx);
        }
        else { // sum_version
          Lambda_m_n(m, N_m(n_idx)) = _sum_version(N_m, lambda_n_m, m, n_idx);
        }
      } // end n_idx
    } // end m

    // Part B: Message passing from symbol nodes to check nodes
    itpp::vec lambda_n(_N);
    for (int n=0; n<_N; ++n) {
      itpp::ivec M_n = _symbol_node_array(n); // set of check nodes connected to symbol node n
      for (int m_idx=0; m_idx<M_n.size(); ++m_idx) {
        _update_lambda_n_m(lambda_n_m,
                           Lambda_m_n, L_u_n, M_n,
                           m_idx,
                           n);
        lambda_n(n) = lambda_n_m(M_n(m_idx), n) + Lambda_m_n(M_n(m_idx), n); // soft bits
      } // end m_idx
    } // end n

    // 3. Make hard decisions
    //
    if (iter > stop_criterion) {
      decoded_bits = itpp::zeros_b(_N);
      for (int n=0; n<_N; ++n) {
        if (lambda_n(n) > 0) {
          decoded_bits(n) = 0;
        }
        else {
          decoded_bits(n) = 1;
        }
      } // end n
      corrected = verify_codeword(decoded_bits); // valid codeword? done?
    }

    iter++;
  } // iteration loop

  decoded_bits.del(_k_b * _Z); // drop parity bits
}


// The product version of belief propagation LDPC decoding
// (as well as the employed linearization) is for example
// described in
//
//    S. Papaharalabos, P. Sweeney, B.G. Evans,
//    P.T. Mathiopoulos, G. Albertazzi, A. Vanelli-Coralli,
//    and G.E. Corazza.
//    "Modified Sum-Product Algorithms for Decoding
//    Low-Density Parity-Check Codes,"
//    IET Communications, vol. 1, no. 3, pp. 294--300,
//    June 2007.
//
const double LDPC_802_11_codec::_product_version(const itpp::ivec & N_m, // neighborhood of m
                                                 const itpp::mat & lambda_n_m,
                                                 const int m,
                                                 const int n_idx) const
{
  double prod = 1.0;
  for (int n_prime_idx=0; n_prime_idx<N_m.size(); ++n_prime_idx) {
    if (n_prime_idx != n_idx) {
      prod *= tanh_lin(0.5 * lambda_n_m(m, N_m(n_prime_idx)));
      // prod *= std::tanh(0.5 * lambda_n_m(m, N_m(n_prime_idx)));  // not as good
      // prod *= tanh_quant(0.5 * lambda_n_m(m, N_m(n_prime_idx))); // not as good
    }
  }
  return 2.0 * atanh_lin(prod); // std::atanh / atanh_quant not as good
  // return 2.0 * std::atanh(prod);  // not as good
  // return 2.0 * atanh_quant(prod); // not as good
}


// The sum version of belief propagation LDPC decoding
// is for example presented in
// 
//    Y.M. Chang, A.I. Vila Casado, M.C. Frank Chang,
//    and R.D. Wesel,
//    "Lower-Complexity Layered Belief-Propagation
//    Decoding of LDPC codes,"
//    in Proc. IEEE International Conference on
//    Communications, Beijing, China, May 2008,
//    pp. 1155--1160.
// 
const double LDPC_802_11_codec::_sum_version(const itpp::ivec & N_m, // neighborhood of m
                                             const itpp::mat & lambda_n_m,
                                             const int m,
                                             const int n_idx) const
{
  double sum = NAN;
  double sign = 1.0;
  for (int n_prime_idx=0; n_prime_idx<N_m.size(); ++n_prime_idx) {
    if (n_prime_idx != n_idx) {
      sum = _soft_xor_exact(sum, std::abs(lambda_n_m(m, N_m(n_prime_idx))));
      // sum = _soft_xor_approx(sum, std::abs(lambda_n_m(m, N_m(n_prime_idx)))); // not as good
      if (lambda_n_m(m, N_m(n_prime_idx)) < 0) {
        sign *= -1;
      }
    }
  }
  return sign * sum;
}


// exact implementation of soft-XOR (see paper
// about sum version listed above)
//
const double LDPC_802_11_codec::_soft_xor_exact(const double x, const double y) const
{
  if (std::isnan(x)) {            // Soft-XOR takes two input arguments.
    return y;                     // By initializing one argument to NaN
  }                               // it can be used for accumulation:
  if (std::isnan(y)) {            // double sum = NaN;
    return x;                     // for (int i=0; i<X; ++i) {
  }                               //    get(y);
  return _phi(_phi(x) + _phi(y)); //    _sum += soft_xor_exact(sum, y);
}                                 // }


// helper fuction used in _soft_xor_exact method above
//
const double LDPC_802_11_codec::_phi(const double x) const
{
  return -std::log(std::tanh(0.5*x)); 
  // return -std::log(tanh_lin(0.5*x)); // not as good (?)
}


// Approximation of soft-XOR presented in 
// 
//    M.M. Mansour and N.R. Shanbhag,
//    "High-Throughput LDPC Decoders,"
//    IEEE Transactions on Very Large Scale Integration
//    Sytems, vol. 11, no. 6, pp. 976--996, Dec. 2003.
//
const double LDPC_802_11_codec::_soft_xor_approx(const double x, const double y) const
{
  if (std::isnan(x)) {            // Soft-XOR takes two args. By initializing
    return y;                     // one to NaN it can be used for accumulation:
  }                               // double sum = NaN;
  if (std::isnan(y)) {            // for (int i=0; i<X; ++i) {
    return x;                     //    get(y);
  }                               //    _sum += soft_xor_approx(sum, y);
  return (std::min(x,y)           // }
  + std::max(0.125 * (5 - 2 * std::abs(x + y)), 0.0)
  - std::max(0.125 * (5 - 2 * std::abs(x - y)), 0.0));
}


// Exact implementation of soft-XOR, including the sign
// of the operands (assuming two operands). This method
// is used in layered belief propagation to simplify the
// implementation.
//
const double LDPC_802_11_codec::_sign_soft_xor_exact(const double x, const double y) const
{
  double sign = 1.0;
  if (x < 0) {
    sign = -1.0;
  }
  if (y < 0) {
    sign *= -1;
  }
  return sign * _phi(_phi(std::abs(x)) + _phi(std::abs(y)));
}


// helper method to compute messages from variable nodes
// check nodes
const void LDPC_802_11_codec::_update_lambda_n_m(itpp::mat & lambda_n_m,
                                                 const itpp::mat & Lambda_m_n,
                                                 const itpp::vec & L_u_n,
                                                 const itpp::ivec & M_n,
                                                 const int m,
                                                 const int n) const
{
  double sum = 0;
  for (int m_prime_idx=0; m_prime_idx<M_n.size(); ++m_prime_idx) {
    if (m_prime_idx != m) {
      sum += Lambda_m_n(M_n(m_prime_idx), n);
    }
  }
  lambda_n_m(M_n(m), n) = L_u_n(n) + sum;
}


// Helper function to find the index of an edge from
// check node m to symbol node n (the index refers to
// the index in the set of edges that emerge from
// check node m). This function is used in the layered
// decoding method below.
//
const int LDPC_802_11_codec::_find_idx(const itpp::ivec & N_m, const int n) const
{
  int idx = 0;
  for (int k=0; k<N_m.size(); ++k) {
    if (N_m(k)==n) {
      idx = k;
      break;
    }
  }
  return idx;
}


// Message passing decoding of an LDPC code based on
// the Zigzag Layered Belief Propagation algorithm from
//
//    Y.M. Chang, A.I. Vila Casado, M.C. Frank Chang,
//    and R.D. Wesel,
//    "Lower-Complexity Layered Belief-Propagation
//    Decoding of LDPC codes,"
//    in Proc. IEEE International Conference on
//    Communications, Beijing, China, May 2008,
//    pp. 1155--1160.
// 
//
void LDPC_802_11_codec::decode_layered(itpp::bvec & decoded_bits,
                                       const itpp::vec & received_symbols,
                                       const int max_iterations,
                                       const int check_iterations)
{
  int stop_criterion = 0;                  // start to _check_
  if (max_iterations > check_iterations) { // if valid codeword
    stop_criterion = check_iterations;     // after stop_criterion
  }                                        // iterations

  // 1. Initialization step:
  //
  itpp::mat lambda_n_m = itpp::zeros(_M, _N);
  itpp::mat Lambda_m_n = itpp::zeros(_M, _N);

  itpp::mat b_m_n = itpp::zeros(_M, _N); // auxiliary variable used in backward recursion
  itpp::mat f_m_n = itpp::zeros(_M, _N); // auxiliary variable used in forward recursion

  itpp::vec L_u_n = _L_c * received_symbols;
  for (int n=0; n<_N; ++n) {
    itpp::vec tmp = L_u_n(n) * itpp::ones(_M);
    lambda_n_m.set_col(n, tmp);
  }

  for (int n=0; n<_N; ++n) {
    itpp::ivec M_n = _symbol_node_array(n); // set of check nodes connected to symbol node n

    for (int m=0; m<M_n.size(); ++m) { // check node idx

      itpp::ivec N_m = _check_node_array(M_n(m)); // set of symbol nodes connected to check node M_n(m)
      int idx = _find_idx(N_m, n);
      if (idx==0) {
        f_m_n(M_n(m), n) = lambda_n_m(M_n(m), n);
      }
      else {
        f_m_n(M_n(m), n) = _sign_soft_xor_exact(f_m_n(M_n(m), N_m(idx-1)),
                                                lambda_n_m(M_n(m), n));
      }
    } // end m
  } // end n

  // 2. Iteration step:
  //
  int iter = 1;
  bool corrected = false;
  while (!(corrected) && (iter<=max_iterations)) {

    if (iter % 2) { // iter is odd
      for (int n=_N-1; n>=0; --n) { // begin backward propagation

        itpp::ivec M_n = _symbol_node_array(n); // set of check nodes connected to symbol node n

        for (int m=0; m<M_n.size(); ++m) { // check node idx
          itpp::ivec N_m = _check_node_array(M_n(m)); // set of symbol nodes connected to check node M_n(m)
          int idx = _find_idx(N_m, n);
          if (idx==0) {
            Lambda_m_n(M_n(m), n) = b_m_n(M_n(m), N_m(idx+1));
          }
          else if (idx==N_m.size()-1) {
            Lambda_m_n(M_n(m), n) = f_m_n(M_n(m), N_m(idx-1));
          }
          else {
            Lambda_m_n(M_n(m), n) = _sign_soft_xor_exact(f_m_n(M_n(m), N_m(idx-1)),
                                                         b_m_n(M_n(m), N_m(idx+1)));
          }
        } // end m

        for (int m=0; m<M_n.size(); ++m) { // check node idx

          _update_lambda_n_m(lambda_n_m,
                             Lambda_m_n, L_u_n, M_n,
                             m,
                             n);

          itpp::ivec N_m = _check_node_array(M_n(m)); // set of symbol nodes connected to check node M_n(m)
          int idx = _find_idx(N_m, n);
          if (idx==N_m.size()-1) {
            b_m_n(M_n(m), n) = lambda_n_m(M_n(m), n);
          }
          else {
            b_m_n(M_n(m), n) = _sign_soft_xor_exact(b_m_n(M_n(m), N_m(idx+1)),
                                                    lambda_n_m(M_n(m), n));
          }
        } // end m
      } // end bwd propagation
    }
    else { // iter is even
      for (int n=0; n<_N; ++n) { // begin forward propagation

        itpp::ivec M_n = _symbol_node_array(n); // set of check nodes connected to symbol node n

        for (int m=0; m<M_n.size(); ++m) { // check node idx
          itpp::ivec N_m = _check_node_array(M_n(m)); // set of symbol nodes connected to check node M_n(m)
          int idx = _find_idx(N_m, n);
          if (idx==0) {
            Lambda_m_n(M_n(m), n) = b_m_n(M_n(m), N_m(idx+1));
          }
          else if (idx==N_m.size()-1) {
            Lambda_m_n(M_n(m), n) = f_m_n(M_n(m), N_m(idx-1));
          }
          else {
            Lambda_m_n(M_n(m), n) = _sign_soft_xor_exact(f_m_n(M_n(m), N_m(idx-1)),
                                                         b_m_n(M_n(m), N_m(idx+1)));
          }
        } // end m

        for (int m=0; m<M_n.size(); ++m) { // check node idx

          _update_lambda_n_m(lambda_n_m,
                             Lambda_m_n, L_u_n, M_n,
                             m,
                             n);

          itpp::ivec N_m = _check_node_array(M_n(m)); // set of symbol nodes connected to check node M_n(m)
          int idx = _find_idx(N_m, n);
          if (idx==0) {
            f_m_n(M_n(m), n) = lambda_n_m(M_n(m), n);
          }
          else {
            f_m_n(M_n(m), n) = _sign_soft_xor_exact(f_m_n(M_n(m), N_m(idx-1)),
                                                    lambda_n_m(M_n(m), n));
          }
        } // end m
      } // end fwd propagation
    }

    // 3. Compute soft bits and make hard decisions
    //
    if (iter > stop_criterion) {
      decoded_bits = itpp::zeros_b(_N);
      itpp::vec lambda_n = L_u_n;
      for (int n=0; n<_N; ++n) {
        itpp::ivec M_n = _symbol_node_array(n); // set of check nodes connected to symbol node n
        double sum = 0;
        for (int m_idx=0; m_idx<M_n.size(); ++m_idx) {
          sum += Lambda_m_n(M_n(m_idx), n);
        }
        lambda_n(n) += sum;
        if (lambda_n(n) > 0) {
          decoded_bits(n) = 0;
        }
        else {
          decoded_bits(n) = 1;
        }
      } // end n
      corrected = verify_codeword(decoded_bits); // valid codeword? done?
    }

    iter++;
  } // iteration loop

  decoded_bits.del(_k_b * _Z); // drop parity bits
}
