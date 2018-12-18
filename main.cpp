#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <thread>
#include <vector>
#include <algorithm>
#include <set>

#include "data_preprosessor.hpp"
#include "scheduler.hpp"
#include "worker.hpp"
#include "pso.hpp"
#include "fft.hpp"
#include "additional.hpp"
#include "barrier.hpp"

#define num_dims 4

// Declarations of static members
Scheduler Worker::scheduler = Scheduler(0);
size_t Worker::num_workers;
size_t Worker::num_workers_left;
size_t Worker::num_iterations;
std::vector<std::vector<std::vector<double> > > Worker::data;
bool PSO::break_simulation = false;

void print_results(std::vector<PSO*> workers);

int main(int argc, char *argv[]) {
	// get input arguments
	std::string file_name(argv[1]);
	std::string log_file_name(argv[2]);
	size_t numofsamples_2 = atoi(argv[3]); // number of samples in frequency domain
	size_t W = atoi(argv[4]); // modes (workers) number
	size_t I = atoi(argv[5]); // iterations number
	size_t P = atoi(argv[6]); // number of particles per swarm
	
	// Use single input file containing freqs, amps, phases obtained using get_data.py script
	std::cout << "Filename: " << file_name << std::endl;
	std::cout << "Log filename: " << log_file_name << std::endl;
	std::cout << "numofsamples_2: " << numofsamples_2 << std::endl;
	std::cout << "Number of workers (swarms): " << W << std::endl;
	std::cout << "Number of iterations: " << I << std::endl; //number of iterations for each worker (swarm)
	std::cout << "Number of particles per swarm: " << P << std::endl;

	// number of samples in time domain
	size_t numofsamples = (numofsamples_2 - 1) * 2;

	// Load data: frequencies, amplitudes, phases
	std::vector<double> freqs(DataStream(file_name, 0, numofsamples_2, 0).get());
	std::vector<double> amps(DataStream(file_name, 1, numofsamples_2, 0).get());
	std::vector<double> phases(DataStream(file_name, 2, numofsamples_2, 0).get());

	std::vector<double> summed(DataStream("/development/mssp18/data/summed_04.csv", 0, numofsamples_2, 0).get());

	// double fs = 5000; //5000 is real for the data
	// requaried sampling frequency is double max expected frequency (i.e. can be defined by number of time samples)
	double fs = numofsamples; //5000 is real for the data
	double tk = numofsamples / fs; // measurement time
	std::cout << "tk: " << tk << std::endl;

	// time values of measurements
	std::vector<double> times(numofsamples);
	get_times(times, numofsamples, 1 / fs);

	// filtered amplitues data for the ranges of interest finding
	std::vector<double> amps_gauss(gaussian_filter(summed, 2, 13));

	// std::vector<size_t> dolinki;
	// findMinimas(amps_gauss, 0, numofsamples_2 - 1, dolinki);

	// for (auto d: dolinki)
	// 	std::cout << d << std::endl;

	// return 0;

	// init min/max values
	std::vector<double> xmin, xmax;
	init_mins(xmin, num_dims);
	init_maxs(xmax, num_dims, amps_gauss, fs, 0, amps_gauss.size());

	Worker::set_params(W, I);
	Worker::init_worker_data(num_dims);

    std::vector<PSO*> workers = std::vector<PSO*>();
	std::vector<std::thread> t_workers(Worker::num_workers);
	for (size_t i = 0; i < Worker::num_workers; ++i) {
		workers.push_back(
				new PSO(P, num_dims, xmin, xmax, &times, &amps, &phases,
						numofsamples, i));
		t_workers[i] = std::thread(&PSO::run, std::ref(*workers[i]));
	}

	// wait for workers to finish searching process
	for (auto &t_w: t_workers)
		t_w.join();

	print_results(workers);

	std::vector<std::vector<double> > xmin_helpers;
	std::vector<std::vector<double> > xmax_helpers;

	std::set<std::pair<size_t, size_t>> myset;
	std::vector<std::vector<double>> founds_filtered;

	for (size_t m = 0; m < Worker::num_workers; ++m) {
		// determine frequency ranges for each of modes
		size_t freq = Worker::data.back()[m][1] + freqs[1]; // TODO: moved all results by one freq sample

		std::vector<size_t> left_middle_right = get_left_middle_right(amps_gauss, numofsamples_2, freq);
		size_t freq_left = left_middle_right[0];
		size_t freq_right = left_middle_right[2];

		// this part removes swarms from already selected regions
		if (myset.insert(std::pair<size_t, size_t>(freq_left, freq_right)).second) {
			founds_filtered.push_back(Worker::data.back()[m]);

			xmin_helpers.push_back(std::vector<double> ({xmin[0], freqs[freq_left], xmin[2], phases[freq_left]}));
			xmax_helpers.push_back(std::vector<double> ());
			init_maxs(xmax_helpers.back(), num_dims, amps_gauss, fs, freq_left, freq_right);
		}
	}

	for (auto w: workers)
		delete w;

	size_t H = founds_filtered.size();
	Worker::set_params(H, I);
	PSO::break_simulation = false;

	// copy data that is associated to specified filtered ranges from main workers to helpers
	// ommit founds that have been discarded in founds_filtered
	Worker::data = std::vector<std::vector<std::vector<double> > >();
	for (size_t m = 0; m < founds_filtered.size(); ++m) {
		Worker::data.push_back(std::vector<std::vector<double> >());
		for (size_t k = 0; k < m+1; ++k)
			Worker::data[m].push_back(founds_filtered[k]);
	}

	std::cout << "Number of helpers: " << Worker::num_workers_left << std::endl;
	std::cout << "Number of iterations: " << Worker::num_iterations << std::endl;

    std::vector<PSO*> helpers = std::vector<PSO*>();
	std::vector<std::thread> t_helpers(H);

	for (size_t i = 0; i < H; ++i) {
		helpers.push_back(
				new PSO(P, num_dims, xmin_helpers[i], xmax_helpers[i], &times, &amps, &phases,
						numofsamples, i, nullptr, true));
		t_helpers[i] = std::thread(&PSO::run, std::ref(*helpers[i]));
	}

	for (auto &t_h: t_helpers)
		t_h.join();

	print_results(helpers);

	prepare_log_file_for_visualisation(log_file_name, H, Worker::data.back(), myset, times, amps, amps_gauss, numofsamples);

	for (auto h: helpers)
		delete h;
	return 0;
}

void print_results(std::vector<PSO*> workers) {
	std::cout << std::endl << "Found parameters" << std::endl;
	std::cout << "----" << std::endl;
	for (size_t i = 0; i < workers.size(); ++i) {
		std::cout << workers[i]->get_id() << " ";
		for (double x: workers[i]->getgbest()) {
			std::cout << x << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "----" << std::endl << std::endl;
}
