#include "pso.hpp"

class Helper: public PSO {
public:

	Helper(	size_t numofparticles,
			size_t numofdims,
			std::vector< std::vector<double> > Xmin,
			std::vector< std::vector<double> > Xmax,
			std::vector<double> *time,
			std::vector<double> *A,
			size_t numofsamples,
			size_t id,
			double c1 = 2, double c2 = 2);
};
