#ifndef SIMULATIONSGP_H
#define SIMULATIONSGP_H
#include "AbstractSimulation.h"

class SimulationSGP : public AbstractSimulation
{
public:
	SimulationSGP(int nStep, double timeStep, std::vector<Eigen::Vector3d> startingPositions, Environment environment,
		std::vector<std::shared_ptr<AbstractCompartment>> compartments = {});

	void generateParticles();
	void saveResults(const std::string &path);
};

#endif