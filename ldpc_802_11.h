#ifndef __ldpc_802_11_codec_h
#define __ldpc_802_11_codec_h

#include <iostream>
#include <iomanip> // only used for debugging with cout << setprecision(X) <<
#include <itpp/itbase.h>

// The following class implements LDPC codes defined by IEEE 802.11,
//
//    http://standards.ieee.org/getieee802/download/802.11-2012.pdf.
//
class LDPC_802_11_codec
{
public:

  //
  // Data types
  //

  enum code_rate_type {
    half,
    two_thirds,
    three_quarters,
    five_sixths
  };

  enum decode_method_type {
    // standard sum-product decoding
    product_version,
    // the product is rewritten as a sum (by going to the log domain)
    sum_version
  };

  LDPC_802_11_codec(const int block_length = 648,
                    const code_rate_type rate = half,
                    const decode_method_type = product_version,
                    const bool debug = false);

  //
  // Member functions
  //

  void encode(itpp::bvec & output_bits,
              const itpp::bvec & input_bits) const;

  // Message passing decoding based on the flooding scheme.
  //
  // decoded_bits: Vector of decision bits.
  //
  // received_symbols: Vector of received symbols.
  //
  // max_iterations: Stop belief propagation (if the loop is still
  //                 going) after max_iterations iterations.
  //
  // check_iterations: Start to check if decoded_bits is a valid
  //                   codeword after check_iterations number of
  //                   iterations, and break iterations if valid
  //                   codeword.
  //
  void decode_flooding(itpp::bvec & decoded_bits,
                       const itpp::vec & received_symbols,
                       const int max_iterations = 40,
                       const int check_iterations = 8);
  
  // Zigzag Layered Belief Propagation (Z-LBP) decoding.
  //
  // decoded_bits: Vector of decision bits.
  //
  // received_symbols: Vector of received symbols.
  // 
  // max_iterations: Stop belief propagation (if the loop is still
  //                 going) after max_iterations iterations.
  // 
  // check_iterations: Start to check if decoded_bits is a valid
  //                   codeword after check_iterations number of
  //                   iterations, and break iterations if valid
  //                   codeword.
  
  void decode_layered(itpp::bvec & decoded_bits,
                      const itpp::vec & received_symbols,
                      const int max_iterations = 40,
                      const int check_iterations = 8);

  // Check if a word of bits is a valid codeword.
  bool verify_codeword(const itpp::bvec & codeword) const;

  int get_block_length() const;

  void get_code_rate(int & k,
                     int & n) const;

  void set_channel_reliability_value(const double noise_variance);

private:

  //
  // Helper member functions
  //
  void _generate_matrix_prototype(const int block_length,
                                  const code_rate_type rate,
                                  const bool print_matrix);

  void _compute_row_x(const bool print_row_x);

  itpp::bvec _compute_lambda_i(const int i,
                               const itpp::bvec m) const;

  void _create_parity_check_matrix(const bool print_matrix);

  void _create_tanner_graph(const bool print_data_structures);

  const double _product_version(const itpp::ivec & neighborhood,
                                const itpp::mat & lambda_n_m,
                                const int m,
                                const int n_idx) const;

  const double _sum_version(const itpp::ivec & neighborhood,
                            const itpp::mat & lambda_n_m,
                            const int m,
                            const int idx) const;

  const void _update_lambda_n_m(itpp::mat & lambda_n_m,
                                const itpp::mat & Lambda_m_n,
                                const itpp::vec & L_u_n,
                                const itpp::ivec & M_n,
                                const int m,
                                const int n) const;

  const double _phi(const double x) const;
  const double _soft_xor_exact(const double x, const double y) const;
  const double _soft_xor_approx(const double x, const double y) const;

  const int _find_idx(const itpp::ivec & N_m, const int n) const;
  const double _sign_soft_xor_exact(const double x, const double y) const;

  //
  // Data types
  //
  itpp::mat _matrix_prototype;
  itpp::bmat _parity_check_matrix;
  itpp::Array<itpp::ivec> _check_node_array;  // a set of check nodes connected to symbol nodes
  itpp::Array<itpp::ivec> _symbol_node_array; // a set of symbol nodes connected to check nodes
  const int _block_length;
  const code_rate_type _rate;
  const decode_method_type _decode_method;
  //
  // The following parameters follow the notation in
  //
  //    Z. Cai, J. Hao, P.H. Tan, S. Sun, and P.S. Chin,
  //    "Efficient encoding of IEEE 802.11n LDPC codes,"
  //    Electronics Letters, vol. 42, no. 25, pp. 1471--1472,
  //    Dec. 2006.
  //
  int _Z;
  int _m_b;
  int _n_b;
  int _k_b;
  int _M; // nrof check nodes
  int _N; // nrof symbol nodes
  int _row_x;
  double _L_c;
};

#endif //__ldpc_802_11_codec_h
