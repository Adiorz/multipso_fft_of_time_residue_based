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
	std::cout << "test" << std::endl;
	// get input arguments
	std::string file_name(argv[1]);

	// Use single input file containing freqs, amps, phases obtained using get_data.py script
	std::cout << "Filename: " << file_name << std::endl;

	std::vector<double> amps(DataStream(file_name, 1, 128, 0).get());
	std::vector<double> phases(DataStream(file_name, 2, 128, 0).get());

	std::vector<double> data;

	ifft(amps, phases, data);

	return 0;
}
