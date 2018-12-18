#pragma once

#include <vector>
#include <mutex>
#include <condition_variable>
#include <chrono>

#include "basic_worker.hpp"
using namespace std::chrono_literals;

class Scheduler {
private:
	size_t num_schedulees;
//      // 5s
//	std::chrono::duration<int64_t, std::milli> timeout = 5 * 1000ms;
        // 0.05s
    std::chrono::duration<int64_t, std::milli> timeout = 5 * 100ms;
	std::vector<std::mutex> ms;
	std::vector<std::condition_variable> cvs;
	std::vector<bool> readys;
	std::vector<bool> processeds;

public:
	Scheduler(size_t num_schedulees);
	Scheduler();

	void send_data_to_worker(BasicWorker *basic_worker, std::vector<double> d);

	void wait_for_worker_to_read(BasicWorker *basic_worker);

	void wait_for_data_and_read(BasicWorker *basic_worker);

};

