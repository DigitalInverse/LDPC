#include <iomanip>
#include <itpp/itcomm.h>
#include "ldpc_802_11.h"

using namespace std;
using namespace itpp;

// Works in the same way as "from:step:to" in MATLAB.
// Apparently not provided by my particular version
// of itpp.
//
template<class T>
Vec<T> linspace_fixed_step(T from, T to, T step = 1)
{
  int points = 0;
  if (0 != step) {
    points = itpp::floor_i(double(to-from)/step)+1;
  }
  if (0 >= points) {
    return Vec<T>(0);
  }

  Vec<T> output(points);
  output(0) = from;
  for (int n = 1; n < points; ++n) {
    output(n) = output(n-1)+step;
  }
  return output;
}

void verify_encoding(const LDPC_802_11_codec & codec)
{
  const int NUMBER = 1000;
  
  //RNG_randomize(); // randomize the random number generator
  
  bool it_works = true;
  for (int ctr = 0; ctr<NUMBER; ++ctr) {
    int k, n;
    codec.get_code_rate(k, n);
    int nrof_info_bits = codec.get_block_length() * k / n;
    bvec input_bits = randb(nrof_info_bits);
    bvec output_bits;
    codec.encode(output_bits, input_bits);
    it_works = codec.verify_codeword(output_bits);
  }
  if (it_works) {
    cout << "efficient encoding algorithm has been implemented correctly" << endl;
  }
  else {
    cout << "efficient encoding algorithm has NOT been implemented correctly " << endl;
  }
}

void simulate_BER(LDPC_802_11_codec & codec,
		  const double SNR_low,
		  const double SNR_high,
		  const double SNR_step = 1.0,
		  const bool flooding = true,
		  const bool benchmark = false,       // the default operation is to evaluate BER
		  const int THRESHOLD = 100,          // stop criterion for Monte Carlo loop
		  const ivec iterations = "10 20 30") // evaluate BER of these iterations
{
  vec EbN0dB = linspace_fixed_step(SNR_low, SNR_high, SNR_step);
  vec EbN0   = inv_dB(EbN0dB);  // calculate Eb/N0 in a linear scale instead of dB
  vec N0     = pow(EbN0, -1.0); // N0 is the variance of the (complex valued) noise
  BPSK bpsk;                    // modulator class
  BERC berc;
  
  int k, n;
  codec.get_code_rate(k, n);
  int block_length = codec.get_block_length();
  int nrof_info_bits = block_length * k / n;
  
  bvec info_bits = zeros_b(nrof_info_bits);
  bvec encoded_bits = zeros_b(block_length);
  vec symbols = zeros(block_length);
  bvec decoded_bits = zeros_b(nrof_info_bits);
  
  // print a header with labels and iterations
  cout << setw(4) << "% SNR/iter\t";
  for (int ctr=0; ctr<iterations.size(); ++ctr) {
    cout << setw(10) << iterations(ctr) << "\t\t";
  }
  cout << "\n" << endl;
  
  for (int EbN0dB_idx=0; EbN0dB_idx<EbN0dB.length(); ++EbN0dB_idx) {
    
    cout << setw(4) << EbN0dB(EbN0dB_idx) << "\t\t";
    
    double noise_var = (1.0 * n / k) * N0(EbN0dB_idx) / 2;
    codec.set_channel_reliability_value(noise_var); // ideal soft scaling
    double noise_std = sqrt(noise_var);
    
    for (int iter_idx=0; iter_idx<iterations.size(); ++iter_idx) {
      
      RNG_reset(0); // same noise realization for each SNR point (smoother BER curve)
      berc.clear();
      
      double ctr = 0;
      while (ctr < THRESHOLD) {
	
	info_bits = randb(nrof_info_bits);
	codec.encode(encoded_bits, info_bits);
	bpsk.modulate_bits(encoded_bits, symbols);
	
	symbols += noise_std * randn(block_length); // AWGN channel
	
	if (flooding) {
	  codec.decode_flooding(decoded_bits, symbols, iterations(iter_idx));
	}
	else { // Zigzag Layered BP
	  codec.decode_layered(decoded_bits, symbols, iterations(iter_idx));
	}
	
	berc.count(info_bits, decoded_bits);
	if (benchmark) {
	  ctr++; // ctr is used for counting blocks
	}
	else {
	  ctr = berc.get_errors(); // ctr is used for counting errors
	}
      } // ctr
      
      double error_rate = berc.get_errorrate();
      
      cout << setprecision(4) << setw(10) << error_rate << "\t\t";
      
    } // iter_idx
    
    cout << endl; // cout << 0.5*erfc(sqrt(EbN0(EbN0dB_idx))) << endl; // uncoded BPSK
    
  } // EbN0dB_idx
}

// compute time difference between to ticks
//
double diffclock(clock_t clock1, clock_t clock2) {
  double diffticks = clock1 - clock2;
  double diff = (diffticks) / double(CLOCKS_PER_SEC);
  return diff;
}

int main()
{
  // BEGIN LIST OF PARAMETERS
  //
  const int block_length = 1944; // half rate code by default
  
  const double SNR_low  = -2.0;
  const double SNR_high =  2.0;
  const double SNR_step =  0.5;
  
  const bool flooding = true; // false => Zigzag Layered BP
  
  const bool benchmark = false; // measure the execution time of the decoding algorithm?
  //
  // END LIST OF PARAMETERS
  
  LDPC_802_11_codec ldpc(block_length,
			 LDPC_802_11_codec::half,
			 LDPC_802_11_codec::product_version); // instantiate codec
  
  //verify_encoding(ldpc);
  
  if (benchmark) {
    clock_t begin = clock();
    simulate_BER(ldpc, SNR_low, SNR_high, SNR_step, flooding, true, 100); // simulate 100 blocks
    clock_t end = clock();
    double exec_time = diffclock(end, begin);
    cout << "\n% Execution time: " << exec_time << "\n" << endl;
  }
  else {
    simulate_BER(ldpc, SNR_low, SNR_high, SNR_step, flooding); // default threshold of 100 errors
    cout << endl;
  }
  
  return 0;
}
