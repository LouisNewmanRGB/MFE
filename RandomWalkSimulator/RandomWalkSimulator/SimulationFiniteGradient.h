#ifndef SIMULATIONFINITEGRADIENT_H
#define SIMULATIONFINITEGRADIENT_H

#include "AbstractSimulation.h"

class SimulationFiniteGradient : public AbstractSimulation
{
private:
	std::vector<double> m_sequenceTimes;
	std::vector<Eigen::Vector3d> m_sequenceValues;

public:
	SimulationFiniteGradient(int nStep, double timeStep, std::vector<Eigen::Vector3d> startingPositions, Environment environment,
		std::vector<std::shared_ptr<AbstractCompartment>> compartments, std::vector<double> sequenceTimes, std::vector<Eigen::Vector3d> sequenceValues);

	void generateParticles();
	void saveResults(const std::string &path);
};

#endif

