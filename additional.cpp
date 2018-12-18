//#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <assert.h>

#include "additional.hpp"
#include "fft.hpp"
#include "pso.hpp"


double sigma(double a, double b) {
	double a_2 = a*a;
	double b_2 = b*b;
	return (a_2 - b_2)/(a_2 + b_2);
}

double get_undamped_amp(double damped_amp, double freq, double damp, double tk) {
	return damped_amp * tk * (2.0 * M_PI * freq * damp) / (1.0 - exp(-2.0 * M_PI * freq * damp * tk));
}

double get_freq_of_max_amp(std::vector<double> &data, size_t fs, size_t f_l, size_t f_r) {
	double max = *std::max_element(data.begin() + f_l, data.begin()  + f_r);
	size_t f_idx = std::distance(data.begin() + f_l, std::find(data.begin() + f_l, data.begin() + f_r, max)) + f_l;

	double f_m = (f_idx + 1) * fs/(data.size() - 1)/2;
	std::cout << "Max amp found at " << f_idx << " (" << f_m << " [Hz])" << std::endl;

	return f_m;
}

void init_maxs(std::vector<double> &xmax, size_t numofdims, std::vector<double> &data, size_t fs, size_t f_l, size_t f_r) {

	xmax = std::vector<double>(numofdims);

	double f_m = get_freq_of_max_amp(data, fs, f_l, f_r);
	double max_damp = 0.3;
	double tk = 2 * data.size() / (double)fs;
	double max_abs = get_undamped_amp(*std::max_element(data.begin() + f_l, data.begin() + f_r), f_m, max_damp, tk);

	xmax[0] = max_abs; // amp
	// xmax[0] = max_abs/10.0; // amp
	// xmax[0] = 100.0; // amp
	xmax[1] = fs / 2; // freq
	xmax[2] = max_damp; // damping
	std::cout << "max_abs: " << xmax[0] << std::endl;
	std::cout << "max_damp: " << max_damp << std::endl;
}

void init_mins(std::vector<double> &xmin, size_t numofdims) {
	xmin = std::vector<double>(numofdims);

	xmin[0] = 0.0; //amp
	xmin[1] = 20; //freq
	xmin[2] = 0.001; //damping	
}

void get_frequencies(std::vector<double> &freq, size_t numofsamples,
		double fs) {
	size_t numofsamples_2 = (size_t) (numofsamples / 2) + 1;

	freq = std::vector<double>(numofsamples_2);

	for (size_t i = 0; i < freq.size(); ++i)
		freq[i] = i * fs / (numofsamples_2 - 1) / 2;
}

void get_times(std::vector<double> &times, size_t numofsamples, double ts) {
	times = std::vector<double>(numofsamples);
	std::generate(times.begin(), times.begin() + numofsamples,
			[ts] () mutable ->double {static double value; return value += ts;});
}

void calc_response(std::vector<std::vector<double>> results,
		size_t numofsamples, double ts, std::vector<double> &response) {
	response = std::vector<double>(numofsamples, 0.);
	for (size_t i = 0; i < results.size(); ++i) {
		double amp = results[i][0];
		double freq = results[i][1];
		double damp = results[i][2];
		for (size_t j = 0; j < numofsamples; ++j) {
			double t = j * ts;
			response[j] += amp
					* sin(2 * M_PI * freq * sqrt(1 - damp * damp) * t)
					* exp(-2 * M_PI * freq * damp * t);
		}
	}
}

void approximate_amp(std::vector<std::vector<double>> factors,
		std::vector<double> &time_real, size_t numofsamples, double *out) {
	std::vector<double> appr(numofsamples, 0.0);
	for (size_t i = 0; i < factors.size(); ++i) {
		double amp = factors[i][0];
		double freq = factors[i][1];
		double damp = factors[i][2];
		for (size_t j = 0; j < numofsamples; ++j) {
			appr[j] += amp * exp(-2 * M_PI * freq * damp * time_real[j])
					* sin(
							2 * M_PI * freq * sqrt(1 - damp * damp)
									* time_real[j]);
		}
	}

	std::vector<double> A_appr;
	std::vector<double> P_appr;
	fft(appr, A_appr, P_appr);
	size_t numofsamples_2 = (size_t) (numofsamples / 2) + 1;
	for (size_t i = 0; i < numofsamples_2; ++i)
		out[i] = A_appr[i];
}

double min(double *values, size_t size) {
	size_t counter = 0;
	double min = values[0];
	while (std::isinf(min)) {
		min = values[counter];
		++counter;
	}
	for (size_t i = counter; i < size; ++i) {
		if (values[i] < min and !std::isinf(values[i]))
			min = values[i];
	}
	return min;
}

void findMinimas(std::vector<double> &v, size_t start, size_t end,
		std::vector<size_t> &idx) {
	idx = std::vector<size_t>();
	for (unsigned int i = start + 1; i < end; ++i) {
		if (v[i] < v[i - 1]) {
			unsigned int j = i;
			while (v[j] == v[j + 1]) {
				++j;
			}
			++j;
			if (v[j] > v[i]) {
				idx.push_back(i);
				i = j;
			}
		}
	}
	if (start == 0 && idx.size() == 0)
		idx.push_back(0);
}

double gaussianDistribution(int x, double mu, double sigma) {
        double d = x - mu;
        double n = 1.0 / (sqrtf(2 * M_PI) * sigma);
        return exp(-d*d/(2 * sigma * sigma)) * n;
}

std::vector<double> gaussian_kernel(double sigma, double mu, size_t size) {
        assert(size%2==1);
        int n = floor(size/2.0);
        std::vector<double> kernel;
        double sum = 0.0;
        for (int i = -n; i < n+1; ++i) {
                double val = gaussianDistribution(i, mu, sigma);
                sum += val;
                kernel.push_back(val);
        }
        for (auto& x: kernel) {
                x /= sum;
        }
        return kernel;
}


std::vector<double> conv(std::vector<double> const &f, std::vector<double> const &g) {
  int const nf = f.size();
  int const ng = g.size();
  int const n  = nf + ng - 1;
  std::vector<double> out(n, double());
  for(auto i(0); i < n; ++i) {
    int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
    int const jmx = (i <  nf - 1)? i            : nf - 1;
    for(auto j(jmn); j <= jmx; ++j) {
      out[i] += (f[j] * g[i - j]);
    }
  }
  return out;
}

std::vector<double> gaussian_filter(std::vector<double> &input, size_t sigma, size_t size) {
	size_t n = size;
	std::vector<double> kernel = gaussian_kernel(sigma, 0.0, n);
	std::vector<double> output = std::vector<double>(input.size());
	size_t cols = input.size();
	size_t kCols = kernel.size();
	size_t kCenterX = kernel.size() / 2;
	for (size_t j = 0; j < cols; ++j)          // columns
			{
		for (size_t n = 0; n < kCols; ++n) // kernel columns
				{
			size_t nn = kCols - 1 - n; // column index of flipped kernel
			// index of input signal, used for checking boundary
			size_t jj = j + (n - kCenterX);

			// ignore input samples which are out of bound
			if (jj >= 0 && jj < cols)
				output[j] += input[jj] * kernel[nn];
		}
	}
	return output;
}

bool should_skip(size_t f, std::vector<std::pair<size_t, size_t>> &to_skip) {
	for (size_t i = 0; i < to_skip.size(); ++i)
		if (f >= to_skip[i].first && f <= to_skip[i].second)
			return true;
	return false;
}

std::string doubleToText(const double & d) {
	std::stringstream ss;
	ss << std::setprecision(std::numeric_limits<double>::max_digits10);
	ss << d;
	return ss.str();
}

std::vector<size_t> get_left_middle_right(std::vector<double>& A_gauss,
		size_t numofsamples_2, size_t f) {
	std::vector<size_t> idxL;
	std::vector<size_t> idxR;
	findMinimas(A_gauss, 0, f, idxL);
	findMinimas(A_gauss, f, numofsamples_2 - 1, idxR);
	if (idxR.size() == 0)
		idxR.push_back(numofsamples_2 - 1);
	return std::vector<size_t> { idxL[idxL.size() - 1] + 1, f, idxR[0] };
}

double get_signal_energy(std::vector<double> signal, size_t start, size_t end) {
	if (start == 1 && end == 0) {
		start = 0;
		end = signal.size();
	}
	//signal energy
	return std::accumulate(signal.begin()+start, signal.begin() + end, 0.0, // start with first element
			[](double x, double y) {return x + y*y;});

}

void prepare_log_file_for_visualisation(std::string log_file_name,
		size_t num_workers, std::vector<std::vector<double> > factors,
		std::set<std::pair<size_t, size_t> > freq_ranges,
		std::vector<double> &time, std::vector<double> &amp,
		std::vector<double> &amp_gauss, size_t numofsamples) {

	size_t numofsamples_2 = size_t(numofsamples / 2) + 1;

	std::ofstream myfile;

	std::vector<std::vector<double>> results;
	std::vector<double> response(numofsamples);

	std::vector<double *> individual_responses(num_workers);
	for (size_t i = 0; i < num_workers; ++i) {
		individual_responses[i] = new double[numofsamples_2];
		approximate_amp(std::vector<std::vector<double> > { factors[i] }, time,
				numofsamples, individual_responses[i]);
	}

	for (size_t i = 0; i < num_workers; ++i) {
		results.push_back(factors[i]);
	}
	for (size_t j = 0; j < numofsamples; ++j) {
		double t = time[j];
		response[j] += PSO::calc_response(results, t);
	}
	std::vector<double> A_found;
	std::vector<double> P;

	fft(response, A_found, P);

	myfile.open(log_file_name, std::ios::out);

	size_t num_dims = factors[0].size();

	myfile << num_workers << " " << num_dims << std::endl;

	std::set<std::pair<size_t, size_t>>::iterator it = freq_ranges.begin();
	for (std::vector<double> factor : factors) {
		for (double x : factor) {
			myfile << doubleToText(x) << " ";
		}
		myfile << doubleToText((*it).first) << " " << doubleToText((*it).second)
				<< " ";
		myfile << std::endl;
		std::advance(it, 1);
	}

	for (size_t i = 0; i < numofsamples_2; ++i) {
		double real_A = amp[i];
		double found = A_found[i];
		myfile << i << ";" << real_A << ";" << found;
		for (size_t h = 0; h < num_workers; ++h)
			myfile << ";" << individual_responses[h][i];
		myfile << ";" << amp_gauss[i];
		myfile << ";" << found-real_A;
		myfile << std::endl;
	}
	myfile.close();

	for (auto r : individual_responses)
		delete[] r;
}

void add_multiplicative_noise(std::vector<double> &input, std::vector<double> &output, double percentage) {
	size_t samples = input.size();

	output = std::vector<double>(samples);

	unsigned int mersenneSeed = 1977;

	std::mt19937_64 generator;
	generator.seed(mersenneSeed);

	std::normal_distribution<double> normal;

	for (int i = 0; i < samples; ++i) {
		output[i] = percentage*input[i]/100.;
	}
	for (int i = 0; i < samples; ++i)
		output[i] = input[i] + output[i]*normal(generator);
}

void add_noise(std::vector<double> &input, std::vector<double> &output, double SNR_dB) {
	size_t samples = input.size();

	output = std::vector<double>(samples);

	// generate initial noise, mean 0.0, stddev 1.0 - standard normal distribution
	const double mean = 0.0;
	const double variance = 1.0;
	std::default_random_engine generator;
	std::normal_distribution<double> dist(mean, variance);
	for (int i = 0; i < input.size(); ++i)
		output[i] = SNR_dB*input[i]/100.;
	for (int i = 0; i < input.size(); ++i)
		output[i] = input[i] + input[i] * output[i]*dist(generator);
}

