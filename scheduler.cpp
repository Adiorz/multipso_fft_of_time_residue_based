#include <chrono>
#include <iostream>

// #include <fstream>

#include "worker.hpp"
#include "scheduler.hpp"

using namespace std::chrono_literals;
using namespace std::chrono;

// std::string file_name = "communication_id2.txt";
// std::ofstream communication_file;

Scheduler::Scheduler(size_t num_schedulees) { //number of processes
	this->num_schedulees = num_schedulees;

    std::chrono::duration<int64_t, std::milli> timeout = this->num_schedulees * 100ms;

	ms = std::vector < std::mutex > (num_schedulees);
	cvs = std::vector < std::condition_variable > (num_schedulees);
	readys = std::vector<bool>(num_schedulees, false);
	processeds = std::vector<bool>(num_schedulees, false);
}

Scheduler::Scheduler(): Scheduler(0) {}

void Scheduler::send_data_to_worker(BasicWorker *basic_worker, std::vector<double> d) {
	size_t i = basic_worker->get_id();
	{
		// if (i == 1) {
		// 	communication_file.open(file_name, std::ios::app);
		// 	nanoseconds ns = duration_cast< nanoseconds >(
		// 		high_resolution_clock::now().time_since_epoch()
		// 	);
		// 	// std::cout << "{\"send_data\": (" << ns.count() << ", ";
		// 	communication_file << "{\"send_data\": (" << ns.count() << ", ";
		// }
		std::lock_guard < std::mutex > lk(ms[i]); //send my data

		//write data to bufer
		((Worker*)basic_worker)->write(d);
		//after writing

		readys[i] = true; //my data
		// if (i == 1) {
		// 	nanoseconds ns = duration_cast< nanoseconds >(
		// 		high_resolution_clock::now().time_since_epoch()
		// 	);
		// 	communication_file << ns.count() << ")}" << std::endl;
		// 	communication_file.close();
		// }
	}
	cvs[i].notify_one(); //notify worker//i+1
}

void Scheduler::wait_for_worker_to_read(BasicWorker *basic_worker) { //i<T-1
	{
		size_t i = basic_worker->get_id();
		// if (i == 1) {
		// 	communication_file.open(file_name, std::ios::app);
		// 	nanoseconds ns = duration_cast< nanoseconds >(
		// 		high_resolution_clock::now().time_since_epoch()
		// 	);
		// 	communication_file << "{\"wait_to_be_read\": (" << ns.count() << ", ";
		// }
		std::unique_lock < std::mutex > lk(ms[i]); // my data
		// if (!cvs[i].wait_for(lk, 30*1000ms, [this, i] {return processeds[i];})) //I (i) am waiting for my data to be read
		if (!cvs[i].wait_for(lk, timeout, [this, i] {return processeds[i];})) //I (i) am waiting for my data to be read
			std::cerr << i << "timed out: processeds[" << i << "] = " << processeds[i] << std::endl;
		// if (i == 1) {
		// 	nanoseconds ns = duration_cast< nanoseconds >(
		// 		high_resolution_clock::now().time_since_epoch()
		// 	);
		// 	communication_file << ns.count() << ")}" << std::endl;
		// 	communication_file.close();
		// }
		processeds[i] = false;
	}
}

void Scheduler::wait_for_data_and_read(BasicWorker *basic_worker) { //for i>0
	size_t i = basic_worker->get_id();
	// if (i == 1) {
	// 	communication_file.open(file_name, std::ios::app);
	// 	nanoseconds ns = duration_cast< nanoseconds >(
	// 		high_resolution_clock::now().time_since_epoch()
	// 	);
	// 	communication_file << "{\"wait_for_data\": (" << ns.count() << ", ";
	// }
	std::unique_lock < std::mutex > lk(ms[i - 1]); // data held by to i-1
	// cvs[i - 1].wait_for(lk, 30*1000ms, [this, i] {return readys[i - 1];}); //I (i) am waiting for data held by to i-1
	cvs[i - 1].wait_for(lk, timeout, [this, i] {return readys[i - 1];}); //I (i) am waiting for data held by to i-1
	readys[i - 1] = false;

	//read data
	((Worker*)basic_worker)->read();
	//after read

	processeds[i - 1] = true; //read data belonging to i-1
	// if (i == 1) {
	// 	nanoseconds ns = duration_cast< nanoseconds >(
	// 		high_resolution_clock::now().time_since_epoch()
	// 	);
	// 	communication_file << ns.count() << ")}" << std::endl;
	// 	communication_file.close();
	// }

	lk.unlock();
	cvs[i - 1].notify_one(); //notify main i-1
}

