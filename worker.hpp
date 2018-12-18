#pragma once

#include <iostream>

#include <vector>
#include "basic_worker.hpp"
#include "scheduler.hpp"

class Worker: public BasicWorker {
	size_t counter = 0;

protected:
	static void write_data(size_t i, std::vector<double> d);
	static void read_data(size_t i);

public:
	static size_t num_workers;
	static size_t num_workers_left;
	static size_t num_iterations;
	static std::vector<std::vector<std::vector<double> > > data;
	static Scheduler scheduler;

	static void set_params(size_t num_of_workers, size_t num_of_iterations) {
		Worker::scheduler = Scheduler(num_of_workers);
		Worker::num_workers = num_of_workers;
		Worker::num_workers_left = num_of_workers;
		Worker::num_iterations = num_of_iterations;
		// data structure: #n - data found by n-th swarm
		// #0
		// #0 #1
		// #0 #1 #2
		// #0 #1 #2 #3
		// etc.
		// e.g.:
		// <0, 0, 0, 0>
		// <0, 0, 0, 0> <0, 0, 0, 0>
		// <0, 0, 0, 0> <0, 0, 0, 0> <0, 0, 0, 0>
		// <0, 0, 0, 0> <0, 0, 0, 0> <0, 0, 0, 0> <0, 0, 0, 0>

		Worker::data = std::vector<std::vector<std::vector<double> > >(num_of_workers);
	}

	static void init_worker_data(size_t num_of_dims) {
		for (size_t i = 0; i < num_workers; ++i) {
			Worker::data[i] = std::vector<std::vector<double> >(i + 1);
			for (size_t j = 0; j < i + 1; ++j) {
				Worker::data[i][j] = std::vector<double>(num_of_dims);
			}
		}
	}

	static void print_data() {
		std::cout << "DATA SIZE: " << data.size() << std::endl;
		for (auto d: data) {
			for (auto dd: d) {
				std::cout << "<";
				for (auto p: dd) {
					std::cout << p << ", ";
				}
				std::cout << "> ";
			}
			std::cout << std::endl;
		}
	}

	Worker(size_t id);
	size_t get_id() const {
		return id;
	}
	void write(std::vector<double>);
	void read();
	void run();
};
